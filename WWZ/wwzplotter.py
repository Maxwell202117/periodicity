#https://github.com/skiehl/wwz
# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""A class for plotting results of the weighted wavelet z-transform analysis.
"""

import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import numpy as np
import os
import sys
import copy
import numpy as np
import matplotlib.gridspec as gridspec
from astropy.timeseries import LombScargle

import matplotlib as mpl

__author__ = "Sebastian Kiehlmann"
__credits__ = ["Sebastian Kiehlmann"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sebastian Kiehlmann"
__email__ = "skiehlmann@mail.de"
__status__ = "Production"
#==============================================================================
# CLASSES
#==============================================================================
TT='J0303-2407'
class WWZPlotter:
    """A class for plotting WWZ results."""
    
    #--------------------------------------------------------------------------
    def __init__(self, wwz, tunit=None):
        """A class for plotting WWZ results.

        Parameters
        ----------
        wwz : wwz.WWZ
            A WWZ instance that is used for plotting.

        Returns
        -------
        None.
        """

        self.wwz = wwz
        self.okay = True

        if wwz.wwz is None:
            print('Note: There is no WWZ transform data stored in this WWZ ' \
                  'instance. There will be nothing to plot.')
            self.okay = False

        if wwz.freq is not None:
            # check if frequencies are linearly scaled:
            freq = self.wwz.freq
            df = np.diff(freq)
            self.linear_freq = np.all(np.isclose(df, df.mean()))
            self.fmin = freq.min()
            self.fmax = freq.max()

            # get periods and check if linearly scaled:
            period = 1. / freq
            dp = np.diff(period)
            self.linear_period = np.all(np.isclose(dp, dp.mean()))
            self.pmin = period.min()
            self.pmax = period.max()

            self.n_ybins = freq.size
        else:
            self.okay = False

        if wwz.tau is not None:
            self.tmin = wwz.tau.min()
            self.tmax = wwz.tau.max()
        else:
            self.okay = False

        if self.okay:
            if self.linear_freq:
                self.ymin = self.fmax
                self.ymax = self.fmin
                self.ymin_alt = self.pmax
                self.ymax_alt = self.pmin
                self.ylabel = f'Frequency [1/{tunit}]' \
                        if isinstance(tunit, str) else 'Frequency'
                self.ylabel_alt = f'Period [{tunit}]' \
                        if isinstance(tunit, str) else 'Period'
            elif self.linear_period:
                self.ymin = self.pmin
                self.ymax = self.pmax
                self.ymin_alt = self.fmin
                self.ymax_alt = self.fmax
                self.ylabel = f'Period [{tunit}]' if isinstance(tunit, str) \
                        else 'Period'
                self.ylabel_alt = f'Frequency [1/{tunit}]' \
                        if isinstance(tunit, str) else 'Frequency'
            else:
                self.ymin = 0
                self.ymax = 1
                self.ymin_alt = 1
                self.ymax_alt = 0
                self.ylabel = 'Non-linear scale'
                self.ylabel_alt = 'Non-linear scale'

    #--------------------------------------------------------------------------
    def _select_map(self, select):
        """Helper method to select a map from a WWZ instance.

        Parameters
        ----------
        select : str
            Select either 'wwz' or 'wwa'.

        Raises
        ------
        ValueError
            Raised if 'select' is not one of the allowed options.

        Returns
        -------
        result : numpy.ndarray
            The selected WWZ or WWA array.

        """

        # check that selection is allowed:
        if select.lower() not in ['wwz', 'wwa']:
            raise ValueError(f"'{select}' is not a valid selection.")

        select = select.lower()
        result = eval(f'self.wwz.{select}')

        # check if result map is available:
        if result is None:
            print(f'No {select.upper()} transform available.')

        result = result.transpose()

        return result

    #--------------------------------------------------------------------------
    def plot_map(
            self, select, ax=None, xlabel='MJD[days]', **kwargs):
        """Plot the resulting map from a WWZ instance.

        Parameters
        ----------
        select : str
            Select either 'wwz' or 'wwa'.
        ax : matplotlib.pyplot.axis, optional
            The axis to plot to. If None is given a new axis is crated. The
            default is None.
        xlabel : str, optional
            The x-axis label. If None is provided no label is placed. The
            default is None.
        kwargs : dict, optional
            Keyword arguments forwarded to the matplotlib.pyplot.imshow()
            function.

        Returns
        -------
        matplotlib.pyplot.axis
            The axis to which the map was plotted.
        matplotlib.image.AxesImage
            The image.
        """

        if not self.okay:
            return None, None

        # select result:
        result = self._select_map(select)
        if result is None:
            return None, None

        # create figure if needed:
        if ax is None:
            __, ax = plt.subplots(1)
        for i in range(0,100):
            result[i,:]=result[i,:]/(np.max(result[i,:]))

       
       # print(np.max(result,axis=1))
        #result[0,:]=result[0,:]/np.max(result[0,:])
        #print(type(result))
        # plot:
        extent = [self.tmin, self.tmax, self.ymin, self.ymax]
        #去除坏点和极端小的点
        result[result<1e-3] = 1e-3
        result[np.isnan(result)]=1e-3
        #print(result)
        im = ax.imshow(
                result, origin='upper', aspect='auto', extent=extent,
                **kwargs)
        import pandas as pd

        data = pd.read_csv('some1.csv')
        from matplotlib.colors import LogNorm
        min_trunc_factor = 100
        ax.contour(data,[0.9544,0.9973], cmap='jet',norm=LogNorm(),extent=extent,origin='lower')
        #sig=np.loadtxt('sig0.95.txt')
        #a=sig[:,0]
        #b=sig[:,1]
        #c=sig[:,2]
        #ax.plot(a,b)
        #plt.scatter(a, b, s=c, color='w', edgecolor='w', alpha=1)
        #plt.colorbar(im,orientation='vertical',fraction=4)
        #print(np.min(result,axis=0))
        # add labels:
        if xlabel:
            ax.set_xlabel(xlabel)
        ax.set_ylabel(self.ylabel)
        #print(np.max(result,axis=1))
        return ax, im

    #--------------------------------------------------------------------------
    def plot_map_avg(
            self, select, statistic='mean', ax=None, ylabel=False, **kwargs):
        """Vertically plot an average along the time axis of the transform map.

        Parameters
        ----------
        select : str
            Select either 'wwz' or 'wwa'.
        statistic : str, optional
            Choose either 'mean' or 'median'. The default is 'mean'.
        ax : matplotlib.pyplot.axis, optional
            The axis to plot to. If None is given a new axis is crated. The
            default is None.
        ylabel : bool, optional
            If True a label is added to the y-axis. The default is False.
        **kwargs : dict
            Keyword arguments forwarded to the matplotlib.pyplot.plot()
            function.

        Raises
        ------
        ValueError
            Raised if 'statistic' is not one of the allowed options.

        Returns
        -------
        matplotlib.pyplot.axis
            The axis to which the data was plotted.
        """

        if not self.okay:
            return None

        # select result:
        result = self._select_map(select)
        if result is None:
            return None, None

        # calculate statistic:
        if statistic not in ['mean', 'median']:
            raise ValueError(f"'{statistic}' is not a valid statistic.")
        elif statistic == 'median':
            result_avg = np.median(result, axis=1)
        else:
            result_avg = np.mean(result, axis=1)

        # create figure if needed:
        if ax is None:
            __, ax = plt.subplots(1)

        # plot:
       # LCCF=np.loadtxt('1y.txt')
        #fig, ax = plt.subplots(dpi=600)
      #  ax.errorbar(LCCF[:,0],LCCF[:,2],fmt='.k',ecolor='k',ms=5,elinewidth=0.5,capsize=0) 
      #  ax.plot(LCCF[:,0],LCCF[:,3],ls=':',color='C1',label='1σ ',linewidth=0.6)
      #  ax.plot(LCCF[:,0],LCCF[:,4],ls='--',color='C2',label='2σ ',linewidth=0.6)
      #  ax.plot(LCCF[:,0],LCCF[:,5],ls='-.',color='C3',label='3σ ',linewidth=0.6)
        


        LCCF=np.loadtxt('t5.txt')
        x1=LCCF[:,0]
        y1=LCCF[:,2]
        y2=LCCF[:,3]
        y3=LCCF[:,4]
        y4=LCCF[:,5]
        #plt.title("y-x scatter")
        #plt.xlabel("x")
        #plt.ylabel("y")
        

        plt.scatter([y for y in y1], x1)
        plt.scatter([y for y in y2], x1,ls=':',color='C1',label='1σ ',s=1)
        plt.scatter([y for y in y3], x1,color='C2',label='2σ ',s=1)
        plt.scatter([y for y in y4], x1,color='C3',label='3σ ',s=1)
        plt.axhline(y=381,label='P=381day',linewidth=2, linestyle='--',color='r')
        plt.axhline(y=602,label='P=602day',linewidth=2, linestyle='--',color='r')
        plt.axhline(y=952,label='P=952days',linewidth=2, linestyle='--',color='r')
        plt.legend()
       # text1="P=641day"
      #  plt.text(120,581,text1)
      #  ax.tick_params(top='on', right='on', which='both')
       # ax.xaxis.set_minor_locator(MultipleLocator(100))
       # ax.yaxis.set_major_locator(MultipleLocator(50))
       # ax.yaxis.set_minor_locator(MultipleLocator(25))
      #  ax.set_ylim(0,200)
      #  ax.set_xlim(0,1000)
    #    ax.set_xlabel('Log[Frequency]')
   #     ax.set_ylabel('Power')
     #   text1="P=641day"
      #  ax.legend()
     #   ax.text(641,140,text1)
        #y = np.linspace(self.ymin, self.ymax, result_avg.size)
        #ax.plot(result_avg[::-1], y, **kwargs)

        # add labels:
        if ylabel:
            ax.set_ylabel(self.ylabel)

        ax.set_xlabel(f'{statistic.capitalize()} {select.upper()}')

        return ax

    #--------------------------------------------------------------------------
    def plot_data(
            self, ax=None, errorbars=True, xlabel=None, ylabel=None, **kwargs):
        """Plot the data stored in a WWZ instance.

        Parameters
        ----------
        ax : matplotlib.pyplot.axis, optional
            The axis to plot to. If None is given a new axis is crated. The
            default is None.
        errorbars : bool, optional
            If True errorbars are shown, if uncertainties were stored in the
            WWZ instance. The default is True.
        xlabel : str, optional
            The x-axis description. If None is provided no label is printed.
            The default is None.
        ylabel : str, optional
            The y-axis description. If None is provided no label is printed.
            The default is None.
        **kwargs : dict
            Keyword arguments forwarded to the matplotlib.pyplot.errorbar()
            function.

        Returns
        -------
        matplotlib.pyplot.axis
            The axis to which the data was plotted.
        """
        #TT='K'
        # check if data is available:
        ls = LombScargle(self.wwz.t, self.wwz.x, self.wwz.s_x)
        phase = (self.wwz.t / 929) % 1
        aa=np.linspace(55800, 58500, 100)
        phase_model = np.linspace(-2, 1, 100)
        best_frequency = 1/929
        mag_model = ls.model(phase_model / best_frequency, best_frequency)
        if self.wwz.t is None:
            print('No data available.')
            return None

        # create figure if needed:
        if ax is None:
            __, ax = plt.subplots(1)

        # plot:
        if errorbars and self.wwz.s_x is not None:
            ax.errorbar(self.wwz.t, self.wwz.x, self.wwz.s_x, **kwargs)
            ax.plot(aa, mag_model, '-k', lw=2,color='b')
        else:
            ax.plot(self.wwz.t, self.wwz.x, **kwargs)
            ax.plot(aa, mag_model, '-k', lw=2)        
        # add labels:
        if isinstance(xlabel, str):
            ax.set_xlabel(xlabel)

        if isinstance(ylabel, str):
            ax.set_ylabel(ylabel)
            ax.set_title(TT,fontsize=25)
        
        
        
        return ax
    #--------------------------------------------------------------------------
         
        
    def plot_colorbar(
            self, ax=None, errorbars=True, xlabel=None, ylabel=None, **kwargs):

        # check if data is available:
        if self.wwz.t is None:
            print('No data available.')
            return None

        # create figure if needed:
        if ax is None:
            __, ax = plt.subplots(1)

        # plot:
        if errorbars and self.wwz.s_x is not None:
            ax.errorbar(self.wwz.t, self.wwz.x, self.wwz.s_x, **kwargs)
        else:
            ax.plot(self.wwz.t, self.wwz.x, **kwargs)

        # add labels:
        if isinstance(xlabel, str):
            ax.set_xlabel(xlabel)

        if isinstance(ylabel, str):
            ax.set_ylabel(ylabel)

       
        
            cmap1 = copy.copy(mpl.cm.viridis)
            norm1 = mpl.colors.Normalize(vmin=0, vmax=100)
            im1 = mpl.cm.ScalarMappable(norm=norm1, cmap=cmap1)
            cbar1 = fig.colorbar(im1, cax=axes[0], orientation='horizontal',ticks=np.linspace(0, 100, 11),label='colorbar with Normalize')
            
        return ax
    
    #--------------------------------------------------------------------------
    def add_right_labels(self, ax):
        """Add ticks and labels to the right side of a plot showing the
        alternative unit, i.e. frequency if period is used on the left side and
        vice versa.

        Parameters
        ----------
        ax : matplotlib.pyplot.axis, optional
            The axis to plot to. If None is given a new axis is crated. The
            default is None.

        Returns
        -------
        ax2 : matplotlib.pyplot.axis
            The new axis to which the labels were added.
        """

        ax2 = ax.twinx()
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()

        ax2.set_ylim(self.ymin_alt, self.ymax_alt)

        sys.stderr = open(os.devnull, "w")  # silence stderr to supress warning
        conversion = lambda x: 1/x
        ax2.set_yscale('function', functions=(conversion, conversion))
        sys.stderr = sys.__stderr__  # unsilence stderr
        ax2.yaxis.set_major_locator(LogLocator(subs='all'))
        ax2.set_ylabel(self.ylabel_alt)

        return ax2

    #--------------------------------------------------------------------------
    def plot(self, select, statistic='mean', errorbars=True, xlabel=None,
             ylabel=None, figsize=None, height_ratios=(2, 1,1),
             width_ratios=(5,2), kwargs_map={}, kwargs_map_avg={},
             kwargs_data={}):
        """Plot the WWZ map, average, and data.

        Parameters
        ----------
        select : str
            Select either 'wwz' or 'wwa'.
        statistic : str, optional
            Choose either 'mean' or 'median'. The default is 'mean'.
        errorbars : bool, optional
            If True errorbars are shown, if uncertainties were stored in the
            WWZ instance. The default is True.
        xlabel : str, optional
            The x-axis description. If None is provided no label is printed.
            The default is None.
        ylabel : str, optional
            The y-axis description. If None is provided no label is printed.
            The default is None.
        figsize : tuple, optional
            Set the figure size. The default is None.
        height_ratios : tuple, optional
            Set the size ratio between the top and bottom panel with two values
            in a tuple. The default is (2, 1).
        width_ratios : tuple, optional
            Set the size ratio between the left and right panel with two values
            in a tuple. The default is (5, 1).
        kwargs_map : dict, optional
            Keyword arguments forwarded to plotting the map. The default is {}.
        kwargs_map_avg : dict, optional
            Keyword arguments forwarded to plotting the map average. The
            default is {}.
        kwargs_data : dict, optional
            Keyword arguments forwarded to plotting the data. The default is
            {}.

        Returns
        -------
        ax_map : matplotlib.pyplot.axis
            The map axis.
        ax_map_avg : matplotlib.pyplot.axis
            The map average axis.
        ax_data : matplotlib.pyplot.axis
            The data axis.
        """
        plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(4,2,hspace=0,wspace=0, height_ratios=[3,5,1,0.5],width_ratios=[2.5,1]) # 将窗口分为三行三列
        ax5=plt.subplot(gs[2, :])
        ax5.set_visible(False)
        ax4 = plt.subplot(gs[3, :])
        ax_map = plt.subplot(gs[1, 0]) # 第一个图占据第一行的所有列
        ax_data = plt.subplot(gs[0, 0]) # 第二个图占据第二行前两列
        ax_map_avg = plt.subplot(gs[1,1]) # 第三个图为第三列从第二行以后的行数
        #ax4 = plt.subplot(gs[-1, 0]) # 第四个图为倒数第一行的第一列
       # ax5 = plt.subplot(gs[-1, -2]) # 第五个图为倒数第一行的倒数第二列
        
 
        
        
        
        
        
         
        
        
        # create figure:
        #plt.figure(figsize=figsize)
       # grid = gs.GridSpec(
         #       3, 2, hspace=0, wspace=0, height_ratios=height_ratios,
        #        width_ratios=width_ratios)
       # ax_map = plt.subplot(grid[0,0])
       # ax_data = plt.subplot(grid[1,0])
      #  ax_map_avg = plt.subplot(grid[0,1])
       # #ax_bar= plt.subplot(grid[2,1])
        #ax_colorbar= plt.subplot(grid[1,1])
        #ax.subplots_adjust(right=0.9)
        #cb=ax.colorbar(im,cax=position)
        #cbar=plt.colobar()

        # plot map:
        self.plot_map(
                select, ax=ax_map, xlabel='MJD[days]', **kwargs_map)
        #plt.colorbar(orientation='horizontal')
        # plot map average:
        self.plot_map_avg(
                select, statistic=statistic, ax=ax_map_avg, **kwargs_map_avg)
        extend = (self.ymax - self.ymin) / (self.n_ybins - 1) / 2.
        ax_map_avg.set_ylim(self.ymin-extend, self.ymax+extend)

         #plot data:
        self.plot_data(
                ax=ax_data, errorbars=errorbars, ylabel=ylabel,
                **kwargs_data)
        #ax.plot(aa, mag_model, '-k', lw=2)
        #ax_data.set_xlim(self.tmin, self.tmax)

       # ax_data.set_xlim(self.tmin, self.tmax)
        cmap = mpl.cm.jet
        norm = mpl.colors.LogNorm(vmin=0.001, vmax=1)

        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),shrink=0.8,
             cax=ax4, orientation='horizontal')
        #cbar.set_label('uniform')
        
        # add right axis labels:
       # self.add_right_labels(ax_map_avg)

        # add data axis labels:
        ax_data.set_xlabel(xlabel)
        ax_data.set_ylabel(ylabel)
        plt.setp(ax_data.get_xticklabels(), visible=False)
        plt.setp(ax_map_avg.get_yticklabels(), visible=False)
        plt.setp(ax_map.get_xticklabels(), visible=True)
        plt.savefig(TT+'.pdf',bbox_inches='tight')
        return  ax_map,ax_map_avg,ax_data,ax4

  #  ax_colorbar
  #  

#==============================================================================
