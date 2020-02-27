# Author: Laura Kulowski


''' format directory to save results in and plot style '''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
import sys, os

#-----------------------------------------------------------------------------------------------
def setup_directory(save_folder):
    
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    
    if not os.path.exists(f'{save_folder}data/'):
        os.makedirs(f'{save_folder}data/')

    if not os.path.exists(f'{save_folder}plots/'):
        os.makedirs(f'{save_folder}plots/')

    if not os.path.exists(f'{save_folder}plots/b_field/'):
        os.makedirs(f'{save_folder}plots/b_field/')

    if not os.path.exists(f'{save_folder}plots/u_phi_construct/'):
        os.makedirs(f'{save_folder}plots/u_phi_construct/')

    if not os.path.exists(f'{save_folder}plots/u_phi/'):
        os.makedirs(f'{save_folder}plots/u_phi/')

    return 


def plt_mode(p_mode = 'paper', f_size = 40):
    
    '''
    specify plot style
    : param p_mode:      options: 'paper' or 'poster'
    :                    white background for 'paper'
    :                    black background for 'poster'
    '''
    
    from matplotlib import rc
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"
    plt.rcParams['mathtext.fontset'] = 'cm'
    fsize = f_size 
    matplotlib.rcParams.update({'font.size': fsize})

    if p_mode == 'poster': 
        plt.rcParams.update({
            "lines.color": "white",
            "patch.edgecolor": "white",
            "text.color": "black",
            "axes.facecolor": "white",
            "axes.edgecolor": "lightgray",
            "axes.labelcolor": "white",
            "xtick.color": "white",
            "ytick.color": "white",
            "grid.color": "lightgray",
            "figure.facecolor": "black",
            "figure.edgecolor": "black",
            "savefig.facecolor": "black",
            "savefig.edgecolor": "black"})
    return 
    
