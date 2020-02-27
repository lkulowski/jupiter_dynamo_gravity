# Author: Laura Kulowski

''' Use Ferraro's law to extend surface zonal flow into the dynamo region '''

import numpy as np
import sys, os
import matplotlib
import matplotlib.pyplot as plt
from itertools import groupby
from operator import itemgetter
import copy
from scipy.interpolate import interp1d

from mods import generate_grid
from mods import plotting

#------------------------------------------------------------------------------------------------------
class Dynamo_Zonal_Flow():
    
    def __init__(self, grid, save_folder):
        
        '''
        : param grid:             grid object
        : param save_folder:      folder to save results in
        : param B_comp:           component of the zonal magnetic field (options: 'Br', 'Bt')
        : param ND:               True/False; include dipolar component or not
        '''
        
        self.save_folder = save_folder
        
        # grid parameters
        self.grid = grid
        self.r_inner = self.grid.r_inner 
        self.r_dynamo = self.grid.r_outer
        self.nt = self.grid.nt
        self.nr = self.grid.nr 

        # B-field at dynamo surface
        self.B_comp = 'Bt'         # theta-component of the magnetic field              
        self.ND = True             # do not include the dipolar component of the magnetic field 

        # Juno zonal gauss coefficients (units T) 
        self.gn0 = np.array([np.nan, 410244.7,  11670.4, 4018.6, -34645.4, -18023.6, -20819.6,
                             598.4, 10059.2, 9671.8, -2299.5]) * (10.**-9)

        # physical units for radius
        self.L = 69911. * 10.**3                 
        self.a = 71492. * 10.**3                             # Jupiter's equatorial radius (m) 
        self.r_dynamo_phys = self.r_dynamo * self.a
        
        # grid and properties
        self.rc = self.grid.cgl_points()
        self.rcs = self.grid.cheby_shift()
        self.theta, self.lats_gl = self.grid.theta_lat_grids()
        self.wr = self.grid.cheby_wr()
        x, self.w = self.grid.gauss_leg_xw() 
        self.p = self.grid.leg_P()
        dp_dz, self.dp_dtheta = self.grid.leg_P_deriv()
        
    def B_component(self):
        
        ''' B_component evaluated at dynamo surface '''
        
        if self.ND:
            arr_n = np.arange(2, 11)
        else:
            arr_n = np.arange(1, 11)

        B_surf = np.zeros(self.nt)
            
        for n in arr_n:
            if self.B_comp == 'Br':
                B_comp_n = (n + 1) * (self.a / self.r_dynamo_phys)**(n + 2) * self.gn0[n] * self.p[n, :]
                B_surf += B_comp_n
                
            if self.B_comp == 'Bt':
                B_comp_n = - (self.a / self.r_dynamo_phys)**(n + 2) * self.gn0[n] * self.dp_dtheta[n, :]
                B_surf += B_comp_n
        
        return B_surf

    def B_component_spec_trans(self, N = 130):
        
        '''
        B_omega = B / sin(theta); spherical harmonic expansion of B_omega 
        : param N:       degree of spherical harmonic expansion  
        '''
        
        B_surf = self.B_component()

        if self.B_comp == 'Bt':
            B_omega = B_surf / np.sin(self.theta)

        if self.B_comp == 'Br': 
            B_omega = B_surf 

        B_omega_ns = np.zeros(N + 1)

        # calculate spectral coefficients
        for n in range(N + 1):
            B_omega_n = 0.5 * (2. * n + 1.) * np.dot(self.w, B_omega * self.p[n, :])
            B_omega_ns[n] = B_omega_n
            
        return B_omega_ns

    def separate_poloidal_field_contours(self, Bp_streams, nlevels):

        '''
        identify four different contour groups: 1. north, 2. south, 3. multiple, 4. free
        : param Bp_streams:            poloidal field streamlines
        : param nlevels:               number of contour levels
        : return g_north, g_south,
        :           g_mult, g_free:    list containing lists of indices for different contours 
        '''
        
        # define contour levels
        arr_levels = np.linspace(Bp_streams.min(), Bp_streams.max(), nlevels)
        
        # identify contour groups 
        g_north = []       # northern contours that touch the inner/outer boundaries once       
        g_south = []       # southern contours that touch the inner/outer boundaries once 
        g_mult = []        # contours with multiple surface connections 
        g_free = []        # contours with no surface connections

        for ii in range(nlevels - 1):
            # select indices in countour ii
            indx_t, indx_rcs = np.where((Bp_streams >= arr_levels[ii]) & (Bp_streams <= arr_levels[ii + 1]))
            V = list(zip(indx_t, indx_rcs))

            # highlight only contour in Bp_stream
            Bp_streams_contour = np.zeros([self.nt, self.nr + 1])
            Bp_streams_contour[indx_t, indx_rcs] = 1.

            contour_groups = self.separate_contours(V, Bp_streams_contour)

            # sort groups
            for cg in contour_groups:
                indx_t_cg, indx_rcs_cg = zip(*cg)

                # free 
                if self.nr not in indx_rcs_cg:
                    g_free.append(cg)

                else:
                    # find indices on surface
                    indx_rcs_cg_array = np.array(list(indx_rcs_cg))
                    indx_theta_cg_array = np.array(list(indx_t_cg))
                    
                    indx_rcs_surf = np.where(indx_rcs_cg_array == self.nr)[0]
                    indx_theta_surf = indx_theta_cg_array[indx_rcs_surf]
                    indx_theta_surf = np.sort(indx_theta_surf) 

                    # find groupings in surface theta
                    indx_theta_groups = []
                    
                    for k, g in groupby(enumerate(indx_theta_surf), lambda x: x[1] - x[0]):
                        indx_theta_groups.append(list(map(itemgetter(1), g)))

                    #print(len(indx_theta_groups))

                    # if there are multiple surface connections, assign contour to g_mult
                    if len(indx_theta_groups) > 1:
                        g_mult.append(cg)

                    # if there is one surface connection, it could be in north, south, or
                    # multiple (i.e., for the patch at the equator)
                    if len(indx_theta_groups) == 1:
                        # the north/south usually connect to the inner r-boundary
                        if 0 in indx_rcs_cg: 
                            theta_contour_mean = np.mean([self.theta[indx_t_ii] for indx_t_ii in indx_t_cg])

                            if theta_contour_mean < np.pi/2.:
                                g_north.append(cg)

                            if theta_contour_mean > np.pi/2.:
                                g_south.append(cg)
                        else:
                            g_mult.append(cg)

        # check that all grid points have been assigned to a group 
        assigned_grid_pts = len(sum(g_north, []) + sum(g_south, []) + sum(g_mult, []) + sum(g_free, []))
        num_grid_pts = (self.nr + 1) * self.nt
        assert assigned_grid_pts == num_grid_pts, 'Unclassified grid points'
        
        return g_north, g_south, g_mult, g_free
    

    def breath_first_search(self, v0, V, A):

        '''
        separate contours using a breath first search        
        '''
        
        visited = []
        stack = [v0]

        while stack:
            v = stack.pop()
            v_x, v_y = v[0], v[1]

            if v not in visited and A[v] == 1:
                visited.append(v)

                if v_x < A.shape[0] - 1:
                    stack.append((v_x + 1, v_y))

                    if v_y < A.shape[1] - 1:
                        stack.append((v_x + 1, v_y + 1))
                    if v_y >= 1:
                        stack.append((v_x + 1, v_y - 1))

                if v_x >= 1:
                    stack.append((v_x - 1, v_y))

                    if v_y < A.shape[1] - 1:
                        stack.append((v_x - 1, v_y + 1))
                    if v_y >= 1:
                        stack.append((v_x - 1, v_y - 1))

                if v_y < A.shape[1] - 1:
                    stack.append((v_x, v_y + 1))

                if v_y >= 1:
                    stack.append((v_x, v_y - 1))

        return visited

    def separate_contours(self, V, A):

        '''
        separate contours using a breath first search
        '''
        
        contour_groups = []
        v0 = V[0]
        unassigned_indices = [v0]

        while unassigned_indices:
            v0 = unassigned_indices[0]

            contour_i = self.breath_first_search(v0, V, A)
            contour_groups.append(contour_i)

            assigned_indices = sum(contour_groups, [])
            unassigned_indices = list(set(V) - set(assigned_indices))
        
        return contour_groups

    def fill_contours(self, g_north, g_south, g_mult, g_free, omega_surf, m_method = 'mean'):

        # initialize omega_phi 
        omega_phi = np.zeros([self.nt, self.nr + 1])
        omega_phi[:] = np.nan

        # order matters when method is north/south
        if m_method == 'north' or m_method == 'mean':
            g_list = [g_north, g_mult, g_south, g_free]
        if m_method == 'south':
            g_list = [g_south, g_mult, g_north, g_free]
            
        # loop over g_list 
        for ii, g_list_ii in enumerate(g_list):

            # loop over contours within contour group 
            for jj, contour_group in enumerate(g_list_ii):

                # contour indices
                indx_t, indx_rcs = zip(*contour_group)
                indx_t = np.array(indx_t)
                indx_rcs = np.array(indx_rcs)

                # find theta values on surface
                indx_rcs_surf = np.where(indx_rcs == self.nr)
                theta_surf = indx_t[indx_rcs_surf]

                if m_method == 'mean':
                    if ii < 3: 
                        omega_surf_vals = omega_surf[theta_surf]
                        omega_phi[indx_t, indx_rcs] = np.mean(omega_surf_vals)
                    
                if m_method == 'north' or m_method == 'south':     
                    # fill in first contour group with mean values 
                    if ii == 0:
                        omega_surf_vals = omega_surf[theta_surf]
                        omega_phi[indx_t, indx_rcs] = np.mean(omega_surf_vals)

                    # fill in multiple connections with northern/southern value 
                    if ii == 1:
                        # find groupings in surface theta
                        surf_theta_double = []
                        theta_surf = np.sort(theta_surf)
                    
                        for k, g in groupby(enumerate(theta_surf), lambda x: x[1] - x[0]):
                            surf_theta_double.append(list(map(itemgetter(1), g)))

                        # fill in the patch
                        if len(surf_theta_double) == 1:
                            omega_surf_vals = omega_surf[theta_surf]
                            omega_phi[indx_t, indx_rcs] = np.mean(omega_surf_vals)

                        # fill in double connections
                        if len(surf_theta_double) == 2:
                            for theta_g in surf_theta_double:
                                mean_theta_g = np.mean(self.theta[theta_g])

                                if (m_method == 'north') & (mean_theta_g < np.pi / 2):
                                    omega_surf_vals = omega_surf[theta_g]
                                    omega_phi[indx_t, indx_rcs] = np.mean(omega_surf_vals)
                                    
                                if (m_method == 'south') & (mean_theta_g > np.pi / 2):
                                    omega_surf_vals = omega_surf[theta_g]
                                    omega_phi[indx_t, indx_rcs] = np.mean(omega_surf_vals)                        

                    # shift omega_surf profile to match the last double connection value 
                    if ii == 2:
                        # copy the current omega_phi
                        if jj == 0:
                            omega_surf_shift = copy.deepcopy(omega_phi[:, -1])
                            indx_nan = np.where(np.isnan(omega_phi[:, -1]))[0]

                            if m_method == 'north':
                                indx_cont = indx_nan[-1] + 1
                                add_const = omega_phi[:, -1][indx_cont] - omega_surf[indx_cont]
                                omega_surf_shift[indx_nan] = omega_surf[indx_nan] + add_const

                            if m_method == 'south':
                                indx_cont = indx_nan[0] - 1
                                add_const = omega_phi[:, -1][indx_cont] - omega_surf[indx_cont]
                                omega_surf_shift[indx_nan] = omega_surf[indx_nan] + add_const                                

                            # check
                            plt.figure(figsize = (8, 10))
                            plt.plot(omega_surf, self.lats_gl, 'k', label = 'omega_surf')
                            plt.plot(omega_phi[:, -1], self.lats_gl, 'b', label = 'omega_phi (adj)')
                            plt.plot(omega_surf_shift, self.lats_gl, 'r--', label = 'omega_phi (shift)')
                            plt.legend(fontsize = 12)
                            plt.xlabel('omega_phi')
                            plt.ylabel('Latitude (Degrees)')
                            plt.savefig(f'{self.save_folder}/plots/u_phi_construct/check_shift.pdf')
                            
                            plt.close() 

                        omega_surf_vals = omega_surf_shift[theta_surf]
                        omega_phi[indx_t, indx_rcs] = np.mean(omega_surf_vals)

                # fill in free contours
                if ii == 3:
                    scl_fac = 1.19
                    
                    # take a theta cross-section at the mean r-value
                    rcs_free_mean = int(indx_rcs.mean())
                    omega_free_cut = omega_phi[:, rcs_free_mean]
                
                    indx_bad = np.isnan(omega_free_cut)
                    indx_good = np.isfinite(omega_free_cut)
    
                    f_free = interp1d(self.lats_gl[indx_good], omega_free_cut[indx_good])
                    free_vals = f_free(self.lats_gl[indx_bad])
                    omega_phi[indx_t, indx_rcs] = scl_fac * np.mean(free_vals)
                                                                                                
        return omega_phi

    def smooth_omega_phi_spec(self, omega_phi_jag, n_cheby, m_leg):
        
        '''
        the omega_phi profile generated by fill_contours is jagged
        we will smooth by representing it up to degree n_cheby in Chebyshev polynomials
        and m_leg in legendre polynomials
        : parma omega_phi_jag:     unsmoothed omega_phi profile
        : n_cheby:                 degree of chebyshev polynomial for smoothing
        : m_leg:                   degree of legendre polynomial for smoothing
        '''

        # smooth in r-direction
        # spectral
        omega_phi_cheby_n = np.zeros([n_cheby + 1, self.nt])

        for tt in range(self.nt): 
            for nn in range(n_cheby + 1):
                f_nn = omega_phi_jag[tt, :] * self.grid.T(nn, self.rc)
                omega_phi_cheby_unnorm = np.dot(self.wr, f_nn)

                if nn == 0 or nn == (n_cheby - 1): 
                    omega_phi_cheby = 1./ np.pi * omega_phi_cheby_unnorm
                else: 
                    omega_phi_cheby = 2. / np.pi * omega_phi_cheby_unnorm

                omega_phi_cheby_n[nn, tt] = omega_phi_cheby

        # physical
        omega_phi_r_smooth = np.zeros([self.nt, self.nr + 1])

        for nn in range(n_cheby + 1):
            omega_phi_r_smooth_nn = np.outer(omega_phi_cheby_n[nn, :], self.grid.T(nn, self.rc))
            omega_phi_r_smooth += omega_phi_r_smooth_nn

        # smooth in theta-direction
        # spectral 
        omega_phi_leg_n = np.zeros([m_leg + 1, self.nr + 1])

        for rr in range(self.nr + 1):
            for nn in range(m_leg + 1):
                f_nn = omega_phi_r_smooth[:, rr] * self.p[nn, :]
                omega_phi_nn = (2. * nn + 1.) / 2. * np.dot(self.w, f_nn)
                omega_phi_leg_n[nn, rr] = omega_phi_nn

        # physical 
        omega_phi_smooth = np.zeros([self.nt, self.nr + 1])
        
        for nn in range(m_leg + 1):
            omega_phi_smooth_nn = np.outer(self.p[nn, :], omega_phi_leg_n[nn, :])
            omega_phi_smooth += omega_phi_smooth_nn

        return omega_phi_smooth
    
    def optimize_smooth_omega_phi(self, omega_phi_jag, omega_surf):

        '''
        find the optimal degree of smoothing so the smoothed angular velocity at the dynamo surface
        best matches the specified angular velocity at the dynamo surface 
        : param omega_phi_jag:      unsmoothed angular velocity profile
        : param omega_surf:         angular velocity profile specified at the dynamo surface   
        '''
        
        n_max = 25
        arr_cheby_deg = []
        arr_leg_deg = []
        arr_L2 = []

        for n_cheby in range(7, n_max + 1):
            for m_leg in range(7, n_max + 1):
                omega_phi_smooth = self.smooth_omega_phi_spec(omega_phi_jag, n_cheby, m_leg)

                # smoothed internal flow 
                omega_phi_surf_smooth = omega_phi_smooth[:, -1]
                u_phi_surf_smooth = omega_phi_surf_smooth * self.r_dynamo * np.sin(self.theta)

                # true flow at surface 
                u_phi_surf = omega_surf * self.r_dynamo * np.sin(self.theta)
                

                # error between smoothed and true surface flows
                L2_error = np.sqrt( sum((u_phi_surf - u_phi_surf_smooth)**2) )

                # record information
                arr_cheby_deg.append(n_cheby)
                arr_leg_deg.append(m_leg)
                arr_L2.append(L2_error)

        # choose n_cheby and m_leg that minimize the L2 error between the smoothed surface flow and
        # the true surface flow
        indx_opt = np.argmin(arr_L2)
        n_cheby_opt = arr_cheby_deg[indx_opt]
        m_leg_opt = arr_leg_deg[indx_opt]
    
        return n_cheby_opt, m_leg_opt

    
    def interior_flow(self, Bp_streams, nlevels = 40, u_rms = 0.1, m_method = 'mean', plt_color = 'r'):
        '''
        calculate the dynamo zonal flow profile   
        : param Bp_streams:        poloidal field streamlines
        : param nlevels:           number of contour levels for the zonal flow
        : param u_rms:             RMS velocity of the flow
        : param m_method:          method used when extending the dynamo surface angular velocity
        :                          along poloidal field streamlines (options: 'mean', 'north', 'south')
        : param plt_color:         plot color
        '''
                
        # set u ~ B_comp
        B_surf_ns = self.B_component_spec_trans()

        if self.B_comp == 'Bt':
            B_surf = np.polynomial.legendre.legval(np.cos(self.theta), B_surf_ns) * np.sin(self.theta)
            u_surf = -B_surf
            omega_surf = u_surf / ( self.r_dynamo * np.sin(self.theta))    
        
        if self.B_comp == 'Bt':
            B_label = r'$B_{\theta}$'
        if (self.B_comp == 'Bt') & (self.ND == True):
            B_label = r'$B_{\theta}^{ND}$'

        if self.B_comp == 'Br':
            B_label = r'$B_{r}$'
        if (self.B_comp == 'Br') & (self.ND == True):
            B_label = r'$B_{r}^{ND}$'

        plt.figure(figsize = (7, 8))
        matplotlib.rcParams.update({'font.size': 15})
        plt.plot(B_surf, self.lats_gl, 'k', label = B_label)
        plt.plot(u_surf, self.lats_gl, color = (0.3, 0.5, 1.), label = r'$u_{\phi}$')
        plt.xlabel('$u_{\phi}$ /' + f'{B_label}')
        plt.ylabel('Latitude (Degrees)')
        plt.yticks(np.linspace(-90., 90., 7))
        plt.ylim([-90., 90.])
        plt.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout() 
        plt.savefig(f'{self.save_folder}/plots/u_phi_construct/u_phi_propto_B_comp.pdf')
        plt.close() 

        # separate poloidal field contours 
        g_north, g_south, g_mult, g_free = self.separate_poloidal_field_contours(Bp_streams, nlevels = nlevels)

        # check contours
        contour_groups_check = [g_north, g_south, g_mult, g_free]
        contour_names = ['North', 'South', 'Multiple', 'Free']
        
        for ii, group_ii in enumerate(contour_groups_check):
            check_matrix = np.zeros([self.nt, self.nr + 1])
            
            for jj, gp in enumerate(group_ii):
                indx_t, indx_rcs = zip(*gp)
                check_matrix[indx_t, indx_rcs] = jj * 2. + 10.

            plotting.Plot_Half_Donut(check_matrix, self.grid).plot_half_donut(ttl = f'Contour Group: {contour_names[ii]}',
                                                                              cb_orientation = 'vertical', fsize = 20)
            plt.savefig(f'{self.save_folder}/plots/u_phi_construct/contour_group_{contour_names[ii].lower()}')
            plt.close()
            
        # extend surface flow along contours
        omega_phi_jag = self.fill_contours(g_north, g_south, g_mult, g_free, omega_surf, m_method = m_method)
        plotting.Plot_Half_Donut(omega_phi_jag, self.grid).plot_half_donut(ttl = '$\omega_{\phi}$ (unsmoothed)',
                                                                           cb_orientation = 'vertical', fsize = 20)
        plt.savefig(f'{self.save_folder}/plots/u_phi_construct/omega_phi_unsmoothed')
        plt.close()

        # smooth interior flow
        n_cheby_opt, m_leg_opt = self.optimize_smooth_omega_phi(omega_phi_jag, omega_surf)
        print(f'cheby = {n_cheby_opt}, leg = {m_leg_opt}')
        omega_phi = self.smooth_omega_phi_spec(omega_phi_jag, n_cheby = n_cheby_opt, m_leg = m_leg_opt)

        #omega_phi = self.smooth_omega_phi_spec(omega_phi_jag, n_cheby = 7, m_leg = 9)    # 14, 9
        plotting.Plot_Half_Donut(omega_phi, self.grid).plot_half_donut(ttl = '$\omega_{\phi}$ (smoothed)',
                                                                       cb_orientation = 'vertical', fsize = 20)
        plt.savefig(f'{self.save_folder}/plots/u_phi_construct/omega_phi_smoothed')
        plt.close()

        
        # u_phi
        rcs_mesh, theta_mesh = np.meshgrid(self.rcs, self.theta)
        u_phi = omega_phi * rcs_mesh * np.sin(theta_mesh)

        # compare original u_phi to modified u_phi
        plt.figure(figsize = (7, 8))
        matplotlib.rcParams.update({'font.size': 15})
        plt.plot(u_surf, self.lats_gl, 'k', label = 'Original')
        plt.plot(u_phi[:, -1], self.lats_gl, 'r--', label = 'Modified')
        plt.xlabel('$u_{\phi}$ (unscaled)')
        plt.ylabel('Latitude (Degrees)')
        plt.yticks(np.linspace(-90., 90., 7))
        plt.ylim([-90., 90.])
        plt.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout() 
        plt.savefig(f'{self.save_folder}/plots/u_phi_construct/modified_u_phi.pdf')
        plt.close()
        
        # scale flow to rms speed
        N = (self.nt) * (self.nr + 1)
        alpha_rms =  1./u_rms * np.sqrt(1./N * np.sum(u_phi**2))
        u_phi_rms = u_phi / alpha_rms

        # compare B_comp to adjusted u_phi
        plt_scl_fac = 1.12 #1.12
        plt.figure(figsize=(9.8 * plt_scl_fac, 11.7 * plt_scl_fac))
        matplotlib.rcParams.update({'font.size': 40})

        lwid = 4
        plt.plot([0., 0.], [-90., 90.], '0.8', linewidth=lwid)
        plt.plot(u_surf / alpha_rms * 100, self.lats_gl, 'k', linewidth = lwid, label = '$u_{\phi} \propto$' + B_label)
        plt.plot(u_phi_rms[:, -1] * 100, self.lats_gl, color = plt_color, linewidth = lwid,
                 label = 'Adjusted $u_{\phi}$', alpha = 0.9)
        plt.xlabel('Velocity (cm s$^{-1}$)')
        plt.ylabel('Latitude (Degrees)')
        plt.xticks([-30., -20., -10., 0., 10., 20., 30.])
        plt.xlim([-33, 33])
        plt.yticks(np.linspace(-90., 90., 7))
        plt.ylim([-90., 90.])
        plt.legend(loc = 'lower right', fontsize = 30)
        plt.subplots_adjust(left = 0.16)
        plt.savefig(f'{self.save_folder}/plots/u_phi/u_surf_adjusted.pdf')

        print(f'RMS velocity = {np.sqrt(np.mean(u_phi_rms**2))}')

        return u_phi_rms

    
