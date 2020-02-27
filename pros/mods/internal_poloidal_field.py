# Author: Laura Kulowski

''' calculate the internal poloidal field that is first and second order continuous with Jupiter's radial field '''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys

#---------------------------------------------------------------------------------------------------------------------
class Internal_Poloidal_Field():
    def __init__(self, r_dynamo, p_form = 'ar2_br3'):

        '''
        : param r_dynamo:    nondimensional radius of the dynamo surface
        : param p_form:      form of the poloidal magnetic field;
        :                    options: 'ar2_br3' for simple polynomial or 'ohmic' for minimal ohmic
        '''
        
        self.r_dynamo = r_dynamo                           
        self.a = 71492. * 10.**3                                       # Jupiter's equatorial radius (m) 
        self.r_dynamo_phys = self.r_dynamo * self.a
        self.L = self.a                                                # nondimensional length scale
        self.p_form = p_form

        # Juno zonal gauss coefficients (units: T)
        self.gn0 = np.array([np.nan, 410244.7,  11670.4, 4018.6, -34645.4, -18023.6, -20819.6, 598.4, 10059.2, 9671.8, -2299.5]) * (10.**-9)

                
    def calc_alpha_beta_n(self):

        ''' calculate alpha_n and beta_n coefficients for spectral poloidal scalar '''

        # initialize solution matrix
        soln = np.zeros([10, 2])

        # solve matrix equation for each degree
        for n in range(1, 11):
            if self.p_form == 'ar2_br3':
                A = np.array([[1., self.r_dynamo_phys],
                              [0., 1.]])
                b = np.array([[1./n * (self.a/self.r_dynamo_phys)**(n + 2.) * self.gn0[n]],
                              [-(n + 2) / n * self.a**(n + 2.) / self.r_dynamo_phys**(n + 3) * self.gn0[n]]])
                
                
            if self.p_form == 'ohmic': 
                A = np.array([[n * (2 * n + 3.) * self.r_dynamo**(n - 1) * self.L**-1, - n * (2 * n + 1.) * self.r_dynamo**(n + 1) * self.L],
                              [n * (2 * n + 3.) * (n - 1) * self.r_dynamo**(n - 2) * self.L**-2, - n * (2 * n + 1.) * (n + 1) * self.r_dynamo**n]])

                b = np.array([[ (self.a / self.r_dynamo_phys)**(n + 2) * self.gn0[n] ],
                              [- (n + 2) * self.a**(n + 2) / (self.r_dynamo_phys**(n + 3)) * self.gn0[n]]])

            A_inv = np.linalg.inv(A)
            x = A_inv @ b

            if self.p_form == 'ohmic':
                x = x / (self.L**n)

            soln[n - 1, :] = x.reshape(-1)

        return soln

    def S_n(self, n, alpha_beta_n, rcs_phys):
        
        '''
        calculate spectral poloidal scalar
        : param n:               spherical harmonic degree 
        : param alpha_beta_n:    alpha_n and beta_n coefficients for degree n
        : param rcs_phys:        dynamo radius in physical units
        '''
        
        if self.p_form == 'ar2_br3':
            s_n = alpha_beta_n[0] * rcs_phys**2 + alpha_beta_n[1] * rcs_phys**3
            
        if self.p_form == 'ohmic':
            s_n = alpha_beta_n[0] * (2. * n + 3.) * rcs_phys**(n + 1.) - alpha_beta_n[1] * (2. * n + 1.) * rcs_phys**(n + 3) 

        return s_n
    
    def poloidal_streams(self, alpha_beta, grid):
        
        '''
        calculate streamlines for poloidal magnetic field
        : param alpha_beta:     alpha_n and beta_n coefficients for spectral poloidal scalar 
        : param grid:           grid object

        '''
        
        dP_dz, dP_dtheta = grid.leg_P_deriv()
        rcs = grid.cheby_shift()
        rcs_phys = rcs * self.a

        # sum over spherical harmonics
        B_p_sum = np.zeros([grid.nt, grid.nr + 1])
        
        for n in range(1, 11):
             B_p_sum_n = - np.outer(dP_dtheta[n, :], self.S_n(n, alpha_beta[n - 1, :], rcs_phys)/rcs_phys) 
             B_p_sum += B_p_sum_n                                                                                                    

        return B_p_sum

    def check_boundary_condition(self, alpha_beta, grid):

        ''' check that boundary conditions at dynamo surface are satisfed '''
        
        matplotlib.rcParams.update({'font.size': 15})

        p = grid.leg_P()
        theta, lats_gl = grid.theta_lat_grids() 

        # initialize internal/external fields + derivatives 
        B_ext = np.zeros(grid.nt)
        dB_ext_dr = np.zeros(grid.nt)
        
        B_int = np.zeros(grid.nt)
        dB_int_dr = np.zeros(grid.nt)
        
        # calculate fields + derivatives 
        for n in range(1, 11):
            # external  
            B_ext_n = (n + 1) * (self.a / self.r_dynamo_phys)**(n + 2) * self.gn0[n] * p[n, :]
            B_ext += B_ext_n
            
            dB_ext_dr_n = - (n + 1) * (n + 2) * self.a**(n + 2) / self.r_dynamo_phys**(n + 3) * self.gn0[n] * p[n, :]
            dB_ext_dr += dB_ext_dr_n

            # internal  
            B_int_n = 1. / (self.r_dynamo_phys**2) * n * (n + 1) * self.S_n(n, alpha_beta[n - 1, :], self.r_dynamo_phys) * p[n, :]
            B_int += B_int_n

            if self.p_form == 'ar2_br3':
                dB_int_dr_n = n * (n + 1) * alpha_beta[n - 1, 1] * p[n, :]
                dB_int_dr += dB_int_dr_n
                
            if self.p_form == 'ohmic':
                dB_int_dr_n = n * (n + 1) * ( alpha_beta[n - 1, 0] * (2. * n + 3) * (n - 1.) * self.r_dynamo_phys**(n - 2) -  alpha_beta[n - 1, 1] * (2. * n + 1) * (n + 1) * self.r_dynamo_phys**n ) * p[n, :]
                dB_int_dr += dB_int_dr_n

        fig, ax = plt.subplots(1, 2, figsize = (10, 6))
        lwid = 3
        ax[0].plot(B_ext, lats_gl, 'k', linewidth = lwid, label = 'External')
        ax[0].plot(B_int, lats_gl, 'r--', linewidth = lwid, label = 'Internal')
        ax[0].set_yticks(np.linspace(-90., 90., 7))
        ax[0].set_ylim([-90., 90.])
        ax[0].set_xlabel(r'$B_{r}$')
        ax[0].set_ylabel('Latitude (Degrees)')
        ax[0].set_title(r'$B_{r}$')

        ax[1].plot(dB_ext_dr, lats_gl, 'k', linewidth = lwid, label = 'External')
        ax[1].plot(dB_int_dr, lats_gl, 'r--', linewidth = lwid, label = 'Internal')
        ax[1].set_yticks(np.linspace(-90., 90., 7))
        ax[1].set_ylim([-90., 90.])
        ax[1].legend()
        ax[1].set_xlabel(r'$\partial B_{r} / \partial r$')
        ax[1].set_ylabel('Latitude (Degrees)')
        ax[1].set_title(r'$\partial B_{r} / \partial r$')

        plt.tight_layout()
        
        return 
        
        
