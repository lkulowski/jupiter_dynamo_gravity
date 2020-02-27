# Author: Laura Kulowski

''' calculate the zonal gravity harmonics produced by a physically motivated model of dynamo zonal flow '''

import numpy as np
import matplotlib.pyplot as plt
from importlib import reload
import pickle
import sys, os


from mods import formatting 
from mods import generate_grid
from mods import internal_poloidal_field
from mods import dynamo_zonal_flow
from mods import gravity_calc
from mods import plotting

# plot formatting
formatting.plt_mode(p_mode = 'paper')

#------------------------------------------------------------------------------------------------------
# create directories to save data/plots
save_folder = f'outputs/minimal_ohmic_dissipation/'
formatting.setup_directory(save_folder)

#------------------------------------------------------------------------------------------------------
# choose poloidal field form
p_form = 'ohmic'  

# create grid
nr = 300
nt = 300 
r_dynamo = 0.9
r_inner = 0.4
r_outer = r_dynamo
grid = generate_grid.Grid(nr, nt, r_inner, r_outer)

#------------------------------------------------------------------------------------------------------
# calculate internal poloidal field 
p_obj = internal_poloidal_field.Internal_Poloidal_Field(r_dynamo, p_form = p_form)
alpha_beta = p_obj.calc_alpha_beta_n()

# plot internal poloidal field 
Bp_streams = p_obj.poloidal_streams(alpha_beta, grid)
plotting.Plot_Half_Donut(Bp_streams, grid).plot_half_donut(v_min = 0., v_max = Bp_streams.max(), contour_field = True,
                                                           cb_orientation = 'vertical', cb_ttl = 'T m', pwr = 4,
                                                           cb_ticks = np.linspace(0, 5., 6) * 10.**4, fmt = '%1.0f')
plt.savefig(f'{save_folder}plots/b_field/poloidal_field.pdf')

# check that boundary conditions at dynamo surface are satisfied
p_obj.check_boundary_condition(alpha_beta, grid)
plt.savefig(f'{save_folder}plots/b_field/check_boundary_cond.pdf')

#------------------------------------------------------------------------------------------------------
# calculate the dynamo zonal flow (RMS velocity: 10 cm/s)
df_obj = dynamo_zonal_flow.Dynamo_Zonal_Flow(grid, save_folder)
u_phi = df_obj.interior_flow(Bp_streams, plt_color = (0.71, 0.13, 0.17))

plotting.Plot_Half_Donut(u_phi * 100, grid).plot_half_donut(v_min = -31, v_max = 25, contour_field = True, nlevels = 20,
                                                            cb_orientation = 'vertical', cb_ttl = 'cm s$^{-1}$',
                                                            cb_ticks = np.linspace(-30., 20., 6) , fmt = "%.0f", fsize = 40)
plt.savefig(f'{save_folder}/plots/u_phi/u_phi.pdf')

#------------------------------------------------------------------------------------------------------
# calculate the gravity harmonics associated with the flow
twe_obj = gravity_calc.Gravity_Calc(u_phi, grid, save_folder = save_folder)
rho_prime = twe_obj.calc_rho_prime()
jns, jns_scl = twe_obj.calc_Jns(rho_prime)

print(jns_scl[[3, 5, 7, 9]] / (10.**-8))

#------------------------------------------------------------------------------------------------------
# save the data
# poloidal field
np.save(f'{save_folder}/data/alpha_beta', alpha_beta)

# zonal flow
np.save(f'{save_folder}/data/u_phi', u_phi)

# density perturbation
np.save(f'{save_folder}/data/rho_prime', rho_prime)

# important grid information
with open(f'{save_folder}/data/grid_class', 'wb') as g:
    pickle.dump(grid, g)

with open(f'{save_folder}/data/df_class', 'wb') as df:
    pickle.dump(df_obj, df)

#------------------------------------------------------------------------------------------------------
# save run parameters
f = open(save_folder + 'output.txt', 'w')
f.write(f'p_form = {p_form} \n')
f.write(f'r_inner = {r_inner} \n')
f.write(f'r_dynamo = {r_dynamo} \n')
f.write(f'nr = {nr} \n')
f.write(f'nt = {nt} \n')
f.close()


plt.close('all')
