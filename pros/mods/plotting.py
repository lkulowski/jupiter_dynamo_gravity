# Author: Laura Kulowski 

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker
from matplotlib.ticker import LogLocator
from matplotlib.table import table


from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)

#-----------------------------------------------------------------------------------------------
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

class Plot_Half_Donut():
    def __init__(self, field, grid):
        self.field = field
        self.grid = grid
                
        self.nr = self.grid.nr
        self.nt = self.grid.nt
        self.r_inner = self.grid.r_inner
        self.r_outer = self.grid.r_outer
        self.rcs = self.grid.cheby_shift()
        theta, lats_gl = self.grid.theta_lat_grids()
        self.theta = theta

    def theta_coord(self):
        '''
        transform theta into plotting coordinates
        '''
        theta_plot = np.zeros_like(self.theta)
    
        for tt in range(self.nt):
            theta_loop = self.theta[tt]
            if (theta_loop >= -0.) and (theta_loop <= np.pi / 2.):
                theta_p = (np.pi / 2.) - theta_loop
                theta_plot[tt] = theta_p
            if (theta_loop > np.pi/2.) and (theta_loop <= np.pi):
                phi_p = theta_loop - (np.pi / 2.)
                theta_p = 2. * np.pi - phi_p 
                theta_plot[tt] = theta_p
        return theta_plot
    
    def setup_axes(self, fig, rect, n_markers):
        """
        With custom locator and formatter.
        Note that the extreme values are swapped.
        """

        tr = PolarAxes.PolarTransform() 

        rs = np.linspace(self.r_inner, self.r_outer, n_markers)  

        ## angle_ticks = [(0, r"$\frac{1}{2}\pi$"),
        ##                (.25*pi, r"$\frac{1}{4}\pi$"),
        ##                (.5*pi, r"$0$"),
        ##                (-0.25*pi, r"$\frac{3}{4}\pi$"),
        ##                (-0.5*pi, r"$\pi$")]

        angle_ticks = []

        grid_locator1 = FixedLocator([v for v, s in angle_ticks])
        tick_formatter1 = DictFormatter(dict(angle_ticks))

        # set the number of r points you want on plot
        grid_locator2 = FixedLocator(rs)
        tick_formatter2 = None

        grid_helper = floating_axes.GridHelperCurveLinear(
            tr, extremes=(0.5 * np.pi, -0.5 * np.pi, self.r_outer, self.r_inner),  
            grid_locator1=grid_locator1,
            grid_locator2=grid_locator2,
            tick_formatter1=tick_formatter1,
            tick_formatter2=tick_formatter2)

        ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)  
        fig.add_subplot(ax1)

        # create a parasite axes whose transData in RA, cz
        aux_ax = ax1.get_aux_axes(tr)

        aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
        ax1.patch.zorder = 0.9  # but this has a side effect that the patch is
        # drawn twice, and possibly over some other
        # artists. So, we decrease the zorder a bit to
        # prevent this.

        return ax1, aux_ax

    def setup_axes2(self, fig, rect, n_markers, lat_min, lat_max, angle_ticks):
        """
        for thin shells, lat subset 
        """

        tr = PolarAxes.PolarTransform() 

        rs = np.linspace(self.r_inner, self.r_outer, n_markers)  

        if not angle_ticks:
            angle_ticks = []

        grid_locator1 = FixedLocator([v for v, s in angle_ticks])
        tick_formatter1 = DictFormatter(dict(angle_ticks))

        # set the number of r points you want on plot
        grid_locator2 = FixedLocator(rs)
        tick_formatter2 = None

        grid_helper = floating_axes.GridHelperCurveLinear(
            tr, extremes=(lat_max/180 * np.pi, lat_min/180 * np.pi, self.r_outer, self.r_inner),     # 70/180
            grid_locator1=grid_locator1,
            grid_locator2=grid_locator2,
            tick_formatter1=tick_formatter1,
            tick_formatter2=tick_formatter2)

        ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
        fig.add_subplot(ax1)

        # format tick labels (see: https://matplotlib.org/mpl_toolkits/axes_grid/users/axisartist.html) 
        ax1.axis["bottom"].major_ticklabels.set_rotation(-90)
        ax1.axis["bottom"].major_ticklabels.set_va("center")
        ax1.axis["bottom"].major_ticklabels.set_pad(12)

        # create a parasite axes whose transData in RA, cz
        aux_ax = ax1.get_aux_axes(tr)

        aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
        ax1.patch.zorder = 0.9  # but this has a side effect that the patch is
        # drawn twice, and possibly over some other
        # artists. So, we decrease the zorder a bit to
        # prevent this.

        return ax1, aux_ax

    def plot_half_donut(self, n_markers = 3, v_min = None, v_max = None, contour_field = False,
                        nlevels = 20, cmap = "RdBu_r", ttl = "", cb_ttl = "", cb_orientation = 'horizontal',
                        num_cb_ticks = None, pwr = None, cb_ticks = [], fmt = "%.1f", fsize = 40):
        
        matplotlib.rcParams.update({'font.size': fsize})
        
        # coordinates for plotting 
        theta_donut = self.theta_coord()
        R1, T1 = np.meshgrid(self.rcs, theta_donut)

        # create figure  
        fig1 = plt.figure(figsize = (10.5, 11.7))
        ax1, aux_ax1 = self.setup_axes(fig1, 111, n_markers)
        v_pad = 0.05

        # min/max values on colorbar 
        if v_min == None and v_max == None:
            v_min = self.field.min()
            v_max = self.field.max()

        # plot data
        if contour_field:
            matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
            levels = np.linspace(v_min, v_max, nlevels)

            cax1 = aux_ax1.contourf(T1, R1, self.field, norm = MidpointNormalize(midpoint=0., vmin = v_min, vmax = v_max),
                                    levels = levels, cmap = cmap, extend = 'both')
            aux_ax1.contour(T1, R1, self.field, colors = 'k', levels = levels, linewidths = 2) 

        else:
            cax1 = aux_ax1.pcolor(T1, R1, self.field, norm = MidpointNormalize(midpoint=0., vmin = v_min, vmax = v_max),
                                  cmap = cmap)

        
        # plot title
        plt.title(ttl)

        # add colorbar
        if not pwr:
            if cb_ticks == []: 
                cb1 = fig1.colorbar(cax1, orientation=cb_orientation, extend = 'both', format = fmt)
            else:
                cb1 = fig1.colorbar(cax1, orientation=cb_orientation, ticks = cb_ticks, extend = 'both', format = fmt)
            cb1.set_label(cb_ttl, rotation=270, labelpad = 45)

        if pwr:
            if cb_ticks == []: 
                cb1 = fig1.colorbar(cax1, orientation=cb_orientation, extend = 'both', format=OOMFormatter(pwr, fformat=fmt, mathText=True))
            else:
                cb1 = fig1.colorbar(cax1, orientation=cb_orientation, ticks = cb_ticks, extend = 'both', format=OOMFormatter(pwr, fformat=fmt, mathText=True))
            cb1.set_label(cb_ttl, rotation=270, labelpad = 45)

                        
        if num_cb_ticks != None:
            tick_locator = ticker.MaxNLocator(nbins=num_cb_ticks)
            cb1.locator = tick_locator
            cb1.update_ticks()
            
        plt.tight_layout()

        return
    
#-----------------------------------------------------------------------------------------------
class Plot_Gravity_Harmonics():
    def __init__(self, f_name, vibes = (1., 0.57, 0.51), label = "N. Equatorial Region"):
        self.f_name = f_name
        self.vibes = vibes
        self.label = label

    def plot_jns(self, table_max_val = 12):
        # load header
        f = open(self.f_name)
        header = f.readline()
        f.close() 
        header = header.split("\t")
        Jn_min = int(header[2][-1])
        Jn_max = Jn_min + len(header) - 3        
        arr_n = np.arange(Jn_min, Jn_max + 1)
        
        # load data 
        data = np.loadtxt(self.f_name, skiprows=1)
        depth_rj = data[0]
        depth_km = data[1]
        data = data[2:]
        
        # make figure 
        fig, ax = plt.subplots(1, 1, figsize = (11, 7))
        ax.set_yscale("log", nonposy = "clip")

        # plot juno data
        box_color = "1"

        arr_n_juno = np.array([3, 5, 7, 9])
        J_ns_juno = np.array([-4.24, -6.89, 12.39, -10.58]) * 10.**-8
        juno_uncert = np.array([0.91, 0.81, 1.68, 4.35]) * 10.**-8

        plt.errorbar(arr_n_juno[J_ns_juno > 0.], J_ns_juno[J_ns_juno > 0], yerr = juno_uncert[J_ns_juno > 0], fmt = '.',
                     color = self.vibes, ecolor='k', elinewidth=3)
        plt.errorbar(arr_n_juno[J_ns_juno < 0], abs(J_ns_juno[J_ns_juno < 0]), yerr = juno_uncert[J_ns_juno < 0], fmt = '.',
                     color = box_color, ecolor='k', elinewidth=3)

        plt.semilogy(arr_n_juno[J_ns_juno > 0], J_ns_juno[J_ns_juno > 0], 'o', markersize = 9, color = 'k', label = 'Juno Value')
        plt.semilogy(arr_n_juno[J_ns_juno < 0], abs(J_ns_juno[J_ns_juno < 0]), 'o', markersize = 9,
                     color = 'k', markerfacecolor=box_color, markeredgewidth = 2, markeredgecolor = 'k')

        # plot data 
        ax.semilogy(arr_n[data > 0.], data[data > 0.], 'o', markersize = 9, color = self.vibes, label = self.label)
        ax.semilogy(arr_n[data < 0.], abs(data[data < 0.]), 'o', markersize = 9, markerfacecolor = "None", markeredgewidth = 2,
                    markeredgecolor = self.vibes)

        # titles/labels
        ax.set_xlabel("Spherical Harmonic Degree")
        ax.set_ylabel(r"$\Delta J_{n}$")
        ax.legend(bbox_to_anchor=(0.98, 1), loc=2, borderaxespad=0., fontsize=17, handletextpad= -0.3, frameon=False)
        ax.set_title('Gravity Harmonics for Geostrophic Flow with \n Penetration Depth $r =$' + f"{depth_rj:.3f}" + "$R_{J}$" + f" ({depth_km:.0f} km)" )


        plt.tight_layout(pad = 2)

        # add table
        #row_labels = [f"$J_{nn}$" + "$(\\times 10^{8})$" for nn in range(Jn_min, table_max_val + 1)]  #header[2:table_max_val + 1]
        row_labels = [r'$\Delta J_{{{0}}}$'.format(nn) + "$(\\times 10^{8})$" for nn in range(Jn_min, table_max_val + 1)]  
        
        table_data = data[0:table_max_val - 2] / (10.**-8)
        table_data = np.round(table_data, 2)
        table_data = ["{:.2f}".format(item) for item in table_data]

        table_vals = np.array(table_data).reshape([-1, 1])
        
        #data = data.reshape([-1, 1])
        #table_vals = data[0:table_max_val - 2, :] / (10.**-8)
        #table_vals = np.round(table_vals, 2)
        
        tb = table(ax, cellText = table_vals, rowLabels = row_labels, colLabels = None, cellLoc = "left", bbox = [1.23, 0.15, 0.15, 0.61])   # x, y, w, h
        cells = tb.get_celld()

        for ii in range(table_vals.shape[0]):
            cells[ii, 0]._loc = "right"
            
        plt.tight_layout(pad = 2)



        return  

    















