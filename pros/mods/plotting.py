# Author: Laura Kulowski

''' make half donut plot ''' 

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

        '''
        make half donut plot
        : param n_markers:         number of radial grid markers
        : param v_min:             minimum value on plot
        : param v_max:             maximum value on plot
        : param contour_field:     True for contour plot, False for colormap
        : param nlevels:           number of contour levels
        : param cmap:              colormap
        : param ttl:               plot title
        : param cb_orientation:    colorbar orientation, 'horizontal' or 'vertical'
        : param num_cb_ticks:      number of ticks on colorbar
        : param pwr:               colorbar in units of x10^pwr
        : param cb_ticks:          list/array of colorbar tick values
        : param fmt:               format number of decimals on colorbar ticks; default is 1 decimal place
        : param fsize:             plot fontsize
        
        '''
        
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



