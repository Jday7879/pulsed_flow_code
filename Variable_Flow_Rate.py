# -*- coding: utf-8 -*-

# Copyright 2020 Jordan Day
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#


# import numpy as np
# import scipy as sp
from scipy import sparse
# from scipy.sparse import linalg
# import matplotlib.pyplot as plt
# import time
from tqdm import tqdm
from Fluid_classes import *
from Plotting_functions import *
import os
from pathlib import Path
import solver
def boundary_init(influent_width):
    global flux, baffle_length, baffle_pairs,Ly, Lx,dx, dy
    psi = np.zeros((nx + 1, ny + 1))
    boundary = np.zeros((nx + 1, ny + 1))  # We set this to 1 on all boundary points
    boundary[0, :] = 1
    boundary[-1, :] = 1
    boundary[:, 0] = 1
    boundary[:, -1] = 1

    psi[0, 0:ny + 3] = flux
    psi[:, -1] = flux

    channel_width = Ly / (baffle_pairs * 2 + 1)
    if channel_width > influent_width:
        in_out_points = influent_width / dy
    else:  # Default influent region to 1 point if too large
        import warnings
        warnings.warn('Influent width was larger than channel width. Influent has been reset to one point.',UserWarning)
        in_out_points = 1

    for bb in np.arange(baffle_pairs):
        boundary[1:round(nx * baffle_length), round(int(2 * bb + 1) * ny * 1 / (2 * baffle_pairs + 1)) - 1] = 1
        psi[1:round(nx * baffle_length), round(int(2 * bb + 1) * ny * 1 / (2 * baffle_pairs + 1)) - 1] = flux
        boundary[round(nx * (1 - baffle_length) + 1):nx,
        round(int(2 * (bb + 1) * ny * 1 / (2 * baffle_pairs + 1))) - 1] = 1

    psi[0, 0:round(int(1 / 2 * ny * 1 / (2 * baffle_pairs + 1) - in_out_points / 2 + 1))] = 0
    psi[-1, ny - round(int(1 / 2 * ny * 1 / (2 * baffle_pairs + 1))):ny + 1] = flux

    for i in np.arange(in_out_points - 1):
        psi[0, round(
            int(1 / 2 * ny * 1 / (2 * baffle_pairs + 1) - in_out_points / 2 + 1 + i) - 2)] = flux / in_out_points * (
                i + 1)
        psi[-1, ny - round(in_out_points / 2) - round(
            int(1 / 2 * ny * 1 / (2 * baffle_pairs + 1) - i))] = flux / in_out_points * (i + 1)

    #boundary[1:-1, 1:-1] *= 0
    #psi[1:-1, 1:-1] *= 0
    return boundary, psi

nx = 300
ny = 300
nz = 1
Lx = 0.32  # length in m
Ly = 0.45  # length in m
Lz = 0.25  # length in m
Lx = 1  # length in m
Ly = 1  # length in m
Lz = 1  # length in m
hrt = 6  # hrt set in any given time units TC is the conversion to seconds
hrt *= TC_hour  # The hydraulic retention time converted into seconds
baffle_length = 91 / 100  # This is the fraction of the tank a baffle takes up in x
baffle_pairs = 6 # Number of baffle pairs (RHS+LHS) = 1 pair.

file_name = 'Y_28_s_100_RT_7_days_HRT_2_baffles_4_diffusion_off' # 'Y_22_s_500_RT_7_days_HRT_12_baffles_1'
# Baffle param

# Runtime for reactor system
RT = 0.5
RT *= TC_day#
dt_max = 8
k = 20  # Store data every k min, 1 is every min, 2 every 2 min
D = 1 / hrt
stream = False#
anode_side = 'all'  # can either be 'all', 'influent', 'effluent'

dx = Lx / nx
dy = Ly / ny
x = np.linspace(0, Lx, nx).T
y = np.linspace(0, Ly, ny).T
[yy, xx] = np.meshgrid(np.linspace(dy / 2, Ly - dy / 2, ny), np.linspace(dx / 2, Lx - dx / 2, nx))

system = domain(Lx, Ly, Lz, nx=nx, ny=ny)

flux = (Lx * Ly) / hrt  # The "area flux" through the system
nxy = nx * ny
nxy_one = (nx + 1) * (ny + 1)
psi = np.zeros((nx + 1, ny + 1))  # This is the stream function
boundary = np.zeros((nx + 1, ny + 1))  # We set this to 1 on all boundary points
boundary[0, :] = 1
boundary[-1, :] = 1
boundary[:, 0] = 1
boundary[:, -1] = 1

edges = boundary
psi[0, 0:ny + 3] = flux
psi[:, -1] = flux
psi, boundary, in_out_points, in_start, out_start = system.influent_effluent_regions(baffle_pairs, baffle_length,
                                                                                     14*dy, psi, boundary, flux)#0.06

bdata = psi[boundary == 1]

file_a = 'hrt' + str(hrt).replace('.', '_') + '_nx' + str(nx) + '_ny' + str(ny) + '_Lx' + str(Lx) + '_Ly' + str(
    Ly) + '_pairs' + str(baffle_pairs) + '_width' + str(np.round(baffle_length, decimals=1)).replace('.',
                                                                                                     '_') + '.csv'
file_x = 'Ux_' + file_a
file_y = 'Uy_' + file_a
data_folder = Path(os.getcwd(), "Output", "Velocity")

#try:
    #print(ux)
    #ux = np.genfromtxt(data_folder / file_x, delimiter=',')
    #uy = np.genfromtxt(data_folder / file_y, delimiter=',')
#except:
    #psi, ux, uy, resid = solver.steady_state(boundary, psi, nx, ny, dx, dy,
    #                                        error=1e-6)  # Using function to determine steady state
    #np.savetxt(data_folder / file_x, ux, fmt='%.18e', delimiter=',')
    #np.savetxt(data_folder / file_y, uy, fmt='%.18e', delimiter=',')
    #ux_i = 0.5 * (ux[0:nx, :] + ux[1:nx + 1, :])
    #uy_i = 0.5 * (uy[:, 0:ny] + uy[:, 1:ny + 1])
    #speed = np.sqrt(ux_i.T ** 2 + uy_i.T ** 2)
    #lw = 10 * speed / speed.max()
    #fig = plt.figure(figsize=(12, 12))
    #plt.streamplot(xx.T, yy.T, ux_i.T, uy_i.T, color=speed, density=1.5, linewidth=1.5, cmap='coolwarm',
#                   arrowstyle='-|>')  # coolwarm')
    #plt.colorbar()
    #plt.xlim((0, Lx))
    #plt.ylim((0, Ly))
    # print('Velocity fields have been determined and saved in text files \n')

# if stream:
#     ux_i = 0.5 * (ux[0:nx, :] + ux[1:nx + 1, :])
#     uy_i = 0.5 * (uy[:, 0:ny] + uy[:, 1:ny + 1])
#     speed = np.sqrt(ux_i.T ** 2 + uy_i.T ** 2)
#     lw = 10 * speed / speed.max()
#     fig, ax = plt.subplots(figsize=(10, 12))
#     strm = ax.streamplot(xx.T, yy.T, ux_i.T, uy_i.T, color=speed, density=4, linewidth=1.5, cmap='coolwarm',
#                          arrowstyle='-|>')
#     cb = plt.colorbar(strm.lines)
#     cb.set_label(label='Velocity (m/s)', fontsize=20)  # size='large')
#     cb.ax.tick_params(labelsize=20)
#     influent_scatter = np.zeros((2, in_out_points))
#     influent_scatter[0, :] = x[1]
#     influent_scatter[1, :] = np.linspace(y[in_start - 1], y[in_start + in_out_points - 1], in_out_points)
#     effluent_scatter = np.zeros((2, in_out_points))
#     effluent_scatter[0, :] = x[-2]
#     effluent_scatter[1, :] = np.linspace(y[out_start - 1], y[out_start + in_out_points - 1], in_out_points)
#     plt.scatter(effluent_scatter[0, :], effluent_scatter[1, :], c='tab:red', s=300 * np.ones((1, 18)),
#                 label='Effluent Region')
#     plt.scatter(influent_scatter[0, :], influent_scatter[1, :], c='tab:green', s=300 * np.ones((1, 18)),
#                 label='Influent Region')
#     plt.legend(loc=2, fontsize=26)
#     plt.title('HRT {} Hours'.format(int(hrt / TC_hour)), fontsize=30)
#     plt.xlabel('x', fontsize=30)
#     plt.ylabel('y', rotation=0, fontsize=30)
#     plt.xticks(fontsize=26)
#     plt.yticks(fontsize=26)
#     plt.xlim((0, Lx))
#     plt.ylim((0, Ly))

#plot_streamlines(xx,yy,ux,uy,title = 'nx = {},ny = {}'.format(int(nx),int(ny)),new_fig = True,stream_density = 8)
#save_figure('Steamlines_nx_{}_ny_{}'.format(int(nx),int(ny)))
plt.close('all')

# %%
external = np.zeros(boundary.shape)
external[0, :] = 1
external[-1, :] = 1
external[:, 0] = 1
external[:, -1] = 1
internal = boundary - external
bio_loc = np.zeros(boundary.shape)
bio_loc[:, 0:ny] += internal[:, 1:ny + 1]
bio_loc[:, 1:ny + 1] += internal[:, 0:ny]
bio_loc = bio_loc[0:nx, 0:ny]
if baffle_length == 0 and baffle_pairs == 0:
    bio_loc[1:-1, -2] = 1
    bio_loc[1:-1, 1] = 1

positional = np.nonzero(np.mean(bio_loc, 0))
switch = np.zeros(bio_loc.shape)
if anode_side == 'influent':
    switch[:, positional[0][0:20:2]] = 1  # :2 to alternate front and back on off 1:20:2 is back, 0:20:2 is front
    bio_loc *= switch
    print('influent facing anodes are active')
elif anode_side == 'effluent':
    switch[:, positional[0][1:20:2]] = 1  # :2 to alternate front and back on off 1:20:2 is back, 0:20:2 is front
    bio_loc *= switch
    print('effluent facing anodes are active')
else:
    anode_side = 'all'
    print('All anodes are active')

bio_number = np.count_nonzero(bio_loc)
# %%

anode_numbers = np.count_nonzero(np.mean(bio_loc, 0))

# Determine anode area based on biofilm and baffle length
Ai = dx * Lz
A = baffle_length * nx * Ai
# A = Lx * Lz #anode area
# Ai = A/nx # area per cell

Vol = Lx * Ly * Lz * 1e3  # Volume in
Voli = dx * dy * Lz * 1e3  # Local volume in L

convert_m2_l = Ai / Voli
convert_l_m2 = Voli / Ai

za = MicrobialPopulation(3000 * np.ones(bio_number)  # np.random.normal(loc = 1000,scale = 10,size = (nx,2)) #initial
                         , 7.9 / TC_day  #7.9 / TC_day  # consumption
                         , 3*1.2 / TC_day  # growth
                         , 0.02  # decay
                         , 20  # sub monod const
                         , 'Anodophilic'  # Defining anodophilic species as class
                         , 5000
                         , mediator_monod=0.2 * 1)
#za.maximum = 5000

Xm = MicrobialPopulation(0 * 0.05 * 512 * np.ones(bio_number)
                         , 8.2 / TC_day
                         , 0.5 / TC_day
                         , 0.02
                         , 80
                         , 'Methanogenic'
                         , 500)

s = Substrate(100 * np.ones((nx, ny)), influent=100, diffusion=1e-9, name='Acetate')
s.current = s.update_influent(baffle_pairs, in_out_points, ny)

m_total = 1  # mg-M / mg-Z Initial mediator

mox = Parameters(0.99 * m_total * np.ones(bio_number), name='Oxidised Mediator')

mred = Parameters(m_total - mox.initial, name='Reduced Mediator')
Ym = 22#.75  # mg-M /mg -s 36#
m = 2  # mol-e / mol-M
gamma = 663400  # mg-M/mol-M
T = 298  # K
j0 = 1e-2  # 1e-2-> almost identical to 1, but much faster run times
BV_full = False

j_min = 1.60
j_max = 1.60  # 1.34#0.64
##############################################
# changed here and j_0 ###
E_min_scaled = j_min * (A * 500)  # /(R*T/(m*F*j0))Anode area times sum of res
E_max_scaled = j_max * (A * 500)  # /(R*T/(m*F*j0)) Anode area times sum of res

#############################################
j_test = 1.3736263736263736#1.5  # .4#2.4#1.6#2.4/10  # 40 j0 = 1e-4
# Full BV stable for hrt = 2, 6, j0 = 1e-4 , J_test = 1.4

Rin = Parameters(0, minimum=7, maximum=7, k=0.006 * A / Vol, name='Internal Resistance')
Rex = 1  # ohm
E_test = 10  # j_test * (A * (Rin.minimum+Rex)) # j_test*(R*T/(m*F*j0)+500*A*(0.92/0.08)/(0.92/0.08))
E_test = j_test*A*(8)
E = Parameters(0, minimum=E_test, maximum=E_test, k=0.0006, name='Voltage')
# E = Parameters(0, minimum=10, maximum=10, k=0.0006, name='Voltage')



pref = gamma / (m * F)
I = Parameters(0, name='Current (Not to be confused with current value)')
s_blank = np.zeros((bio_number, 100000))
s_blank1 = np.zeros((100000))

ij = 0
t = GeneralVariable(0, name='Time')
Rin.current = Rin.minimum  # +(Rin.maximum - Rin.minimum)*np.exp(-Rin.k*sum(Za.initial)/nx) + Rex
Rin.storage[0] = Rin.current

# setting up locations for biofilm
# bio_loc = np.zeros((nx,ny))
bio_upper = np.zeros((nx, ny))
bio_lower = np.zeros((nx, ny))
bio_lower[:, -2] = 1  # used for plotting
bio_upper[:, 1] = 1
consump = np.zeros((nx, ny))  # setting consumption array
med_dist = np.zeros(consump.shape)

#ux_max = np.max(ux)  # max velocity based on steady state
#uy_max = np.max(uy)  # max vel from steady state
# Creating sparse matrix for biofiolm diffusion]
positions = [-1, 0, 1]
diag_x = np.array([[1 / (dx ** 2)], [-2 / (dx ** 2)], [1 / (dx ** 2)]]).repeat(nx, axis=1)
diag_y = np.array([[1 / (dy ** 2)], [-2 / (dy ** 2)], [1 / (dy ** 2)]]).repeat(ny, axis=1)
Dx = sp.sparse.spdiags(diag_x, positions, nx, nx)  # d/dx mat Alternate approach to using diffuse_S array
Dy = sp.sparse.spdiags(diag_y, positions, ny, ny)  # d/dy mat
kappa_bio = 1e-12  # diffusion rate for biofilm
Dx_bio = sp.sparse.spdiags(diag_x, positions, nx, nx).tolil()
Dx_bio[0, -1] += 1 / (dx ** 2)  # Periodic Boundary Conditions
Dx_bio[-1, 0] += 1 / (dx ** 2)  # Periodic BC
Dx_bio *= kappa_bio  # setting up diffusion array for biofilm
Dx_bio = Dx_bio.tocsr()

Bx = sp.sparse.spdiags(diag_x, positions, nx, nx)  # tolil()
By = sp.sparse.spdiags(diag_y, positions, ny, ny)
Iy = sp.sparse.eye(ny)
Ix = sp.sparse.eye(nx)
Diffuse_s = (sp.sparse.kron(Iy, Bx) + sp.sparse.kron(By, Ix)).tolil()

bio_diffusion_x = sp.sparse.kron(Iy, Bx).tolil()
temp_location = np.zeros((nx + 1, ny + 1))
temp_location[:-1, :-1] = bio_loc

for ii in np.arange(nxy):
    ix = int(ii % nx)
    iy = int(np.floor(ii / nx))
    jj = iy * (nx + 1) + ix
    if boundary[ix, iy] * boundary[ix, iy + 1] == 1:  # Boundary on left
        Diffuse_s[ii, ii] = Diffuse_s[ii, ii] + 1 / (dx ** 2)  #
        if ix != 0:
            Diffuse_s[ii, ii - 1] = 0
    if boundary[ix + 1, iy] * boundary[ix + 1, iy + 1] == 1:  # Boundary on right
        Diffuse_s[ii, ii] = Diffuse_s[ii, ii] + 1 / (dx ** 2)
        if ix != nx - 1:
            Diffuse_s[ii, ii + 1] = 0

    if boundary[ix, iy] * boundary[ix + 1, iy] == 1:  # Boundary below
        Diffuse_s[ii, ii] = Diffuse_s[ii, ii] + 1 / (dy ** 2)
        if iy != 0:
            Diffuse_s[ii, ii - nx] = 0

    if boundary[ix, iy + 1] * boundary[ix + 1, iy + 1] == 1:  # Boundary above
        Diffuse_s[ii, ii] = Diffuse_s[ii, ii] + 1 / (dy ** 2)
        if iy != ny - 1:
            Diffuse_s[ii, ii + nx] = 0
    if temp_location[ix, iy] * temp_location[ix + 1, iy] == 0:
        bio_diffusion_x[ii, ii] += 1 / (dx ** 2)
    if ix == 0:
        bio_diffusion_x[ii, ii] += 1 / (dx ** 2)
    if ix != 0:
        if temp_location[ix - 1, iy] * temp_location[ix, iy] == 0:
            bio_diffusion_x[ii, ii] += 1 / (dx ** 2)

    if temp_location[ix, iy] + temp_location[ix + 1, iy] == 0:
        bio_diffusion_x[ii, ii] -= 1 / (dx ** 2)

Diffuse_s = Diffuse_s.tocsr()
LU = sp.sparse.linalg

bio_diffusion_x = 0*1e-6 * bio_diffusion_x.tocsr()

za.calculate_positional_distribution(bio_loc)
temp = bio_diffusion_x.dot(np.reshape(za.positional_distribution.T, nxy))
temp = np.reshape(temp.T, (ny, nx)).T
temp[bio_loc != 1] = 0  # Deals with mass produced outside of biofilm! (temp fix)

del temp_location

s.storage[:, :, 0] = s.initial
s.current[:, :] = s.initial  # This line causes s.initial to be linked to s.now

#dt = min(dt_max, 1 / (ux_max / dx + uy_max / dy), (dx ** 2 * dy ** 2) / (
#        2 * s.diffusion * (dx ** 2 + dy ** 2)))
# dt = np.floor(dt*100)/100

ii = 0
rk = np.array([[0, 1 / 2, 1 / 2, 1], [0, 1, 1, 1], [1, 2, 2, 1]]).T

# These have been taken out of loop as they remain constant and independent of Z

bound_out = np.zeros(s.current.shape)
bound_in = np.zeros(s.current.shape)
bound_out[-1, :] = 1
bound_out[:, -1] = 1
bound_out[:, 0] = 1
bound_in[-2, :] = 1
bound_in[:, -2] = 1
bound_in[:, 1] = 1

E.current = E.minimum
Rin.current = Rin.minimum
Rsig = Rin.current + Rex

mred.current = m_total - mox.current
eta_conc = R * T / (m * F) * np.log(m_total / mred.current)

med_dist[bio_loc == 1] = Ai * mred.current / mox.current
summation = np.array([Rsig * np.sum(med_dist, 0), ] * nx)
summation_shaped = summation[bio_loc == 1]
j = (mred.current / mox.current * (E.current - eta_conc)) / (
        R * T / (m * F * j0) + summation_shaped)
del eta_conc

print("System will simulate a {} baffle system with a fluid HRT of {} hours and bio runtime of {} days".format(
    2 * baffle_pairs, hrt / TC_hour, RT / TC_day))

start_time_bio = time.time()
storage_steps = 30#int(k * 20)

s.current = s.update_influent(baffle_pairs, in_out_points, ny)
#total_mass_in = (ux[0, (in_start - 1):in_start + in_out_points - 1] * s.initial[0, (in_start - 1):in_start + in_out_points - 1] * Voli / dx).sum()*dt
total_mass_out = 0
total_removed = 0
total_remaining = 0
total_time = time.time()
pbar = tqdm(total=101, desc="Progress of simulation", ncols=100, )

#flux = psi[-1,-1]
flux_maximum = 1*flux
start_time = time.time()
nxy = nx * ny
nxy_one = (nx + 1) * (ny + 1)

nu = 1e-6  # fluid kinematic viscosity
dt_max_diffusion = (dx ** 2 * dy ** 2) / (2 * nu * (dx ** 2 + dy ** 2))
vorticity = np.zeros((nx, ny))  # This is vorticity
rhs = np.zeros((nx + 1, ny + 1))  # This will be the Laplacian of psi
fx = np.zeros((nx + 1, ny))  # Rightward vorticity flux
fy = np.zeros((nx, ny + 1))  # Upward vorticity flux
d_vorticity_dt1 = np.zeros((nx, ny))  # initalised for rk4 algorithm
d_vorticity_dt2 = np.zeros((nx, ny))
diagonal_x = np.array([[1 / (dx ** 2)], [-2 / (dx ** 2)], [1 / (dx ** 2)]]).repeat(nx + 1, axis=1)
diagonal_y = np.array([[1 / (dy ** 2)], [-2 / (dy ** 2)], [1 / (dy ** 2)]]).repeat(ny + 1, axis=1)
positions = [-1, 0, 1]

Ax = sp.sparse.spdiags(diagonal_x, positions, nx + 1, nx + 1)
Ay = sp.sparse.spdiags(diagonal_y, positions, ny + 1, ny + 1)
Iy = sp.sparse.eye(ny + 1)
Ix = sp.sparse.eye(nx + 1)
AA = sp.sparse.kron(Iy, Ax) + sp.sparse.kron(Ay, Ix)
bd = boundary.T
bd = bd.reshape((nxy_one))
AA = AA.tolil()  # converting to alternate form of sparse matrix
st_loop = time.time()
for i in np.arange(nxy_one):
    if bd[i] == 1:
        AA[i, :] = 0
        AA[i, i] = 1
AA = AA.tocsr()
LU = sp.sparse.linalg.factorized(AA.tocsc())
print('AA matrix defined', time.time() - st_loop)

diagonal_x = np.array([[1 / (dx ** 2)], [-2 / (dx ** 2)], [1 / (dx ** 2)]]).repeat(nx, axis=1)
diagonal_y = np.array([[1 / (dy ** 2)], [-2 / (dy ** 2)], [1 / (dy ** 2)]]).repeat(ny, axis=1)
Bx = sp.sparse.spdiags(diagonal_x, positions, nx, nx)
By = sp.sparse.spdiags(diagonal_y, positions, ny, ny)
Iy = sp.sparse.eye(ny)
Ix = sp.sparse.eye(nx)
BB = sp.sparse.kron(Iy, Bx) + sp.sparse.kron(By, Ix)
CC = sp.sparse.csr_matrix((nxy, nxy_one))
CC = CC.tolil()  # Faster to use this form for indexing vs csr format
BB = BB.tolil()

for ii in np.arange(nxy):
    ix = int((ii) % nx)
    iy = int(np.floor((ii) / nx))
    jj = (iy) * (nx + 1) + ix
    if boundary[ix, iy] * boundary[ix, iy + 1] == 1:  # Boundary on left
        BB[ii, ii] = BB[ii, ii] - 2 / (dx ** 2)
        BB[ii, ii + 1] = BB[ii, ii + 1] * 4 / 3
        if ix != 1:
            BB[ii, ii - 1] = 0
        CC[ii, jj] = CC[ii, jj] + 8 / 3 / (dx ** 4)
        CC[ii, jj + nx + 1] = CC[ii, jj + nx + 1] + 8 / 3 / (dx ** 4)
        CC[ii, jj + 1] = CC[ii, jj + 1] - 8 / 3 / (dx ** 4)
        CC[ii, jj + nx + 2] = CC[ii, jj + nx + 2] - 8 / 3 / (dx ** 4)

    if boundary[ix + 1, iy] * boundary[ix + 1, iy + 1] == 1:  # Boundary on right
        BB[ii, ii] = BB[ii, ii] - 2 / (dx ** 2)
        BB[ii, ii - 1] = BB[ii, ii - 1] * 4 / 3
        if ix != nx - 1:
            BB[ii, ii + 1] = 0
        CC[ii, jj] = CC[ii, jj] - 8 / 3 / (dx ** 4)
        CC[ii, jj + nx + 1] = CC[ii, jj + nx + 1] - 8 / 3 / (dx ** 4)
        CC[ii, jj + 1] = CC[ii, jj + 1] + 8 / 3 / (dx ** 4)
        CC[ii, jj + nx + 2] = CC[ii, jj + nx + 2] + 8 / 3 / (dx ** 4)

    if boundary[ix, iy] * boundary[ix + 1, iy] == 1:  # Boundary below
        BB[ii, ii] = BB[ii, ii] - 2 / (dy ** 2)
        BB[ii, ii + nx] = BB[ii, ii + nx] * 4 / 3
        if iy != 1:
            BB[ii, ii - nx] = 0
        CC[ii, jj] = CC[ii, jj] + 8 / 3 / (dy ** 4)
        CC[ii, jj + nx + 1] = CC[ii, jj + nx + 1] - 8 / 3 / (dy ** 4)
        CC[ii, jj + 1] = CC[ii, jj + 1] + 8 / 3 / (dy ** 4)
        CC[ii, jj + nx + 2] = CC[ii, jj + nx + 2] - 8 / 3 / (dy ** 4)

    if boundary[ix, iy + 1] * boundary[ix + 1, iy + 1] == 1:  # Boundary above
        BB[ii, ii] = BB[ii, ii] - 2 / (dy ** 2)
        BB[ii, ii - nx] = BB[ii, ii - nx] * 4 / 3
        if iy != ny - 1:
            BB[ii, ii + nx] = 0
        CC[ii, jj] = CC[ii, jj] - 8 / 3 / (dy ** 4)
        CC[ii, jj + nx + 1] = CC[ii, jj + nx + 1] + 8 / 3 / (dy ** 4)
        CC[ii, jj + 1] = CC[ii, jj + 1] - 8 / 3 / (dy ** 4)
        CC[ii, jj + nx + 2] = CC[ii, jj + nx + 2] + 8 / 3 / (dy ** 4)
print('BB and CC matrix defined')
CC = CC.tocsr()
BB = BB.tocsr()
del AA, Ax, Bx, Ay, By, diagonal_x, diagonal_y, Ix, Iy

bdata = psi[boundary == 1]
rhs[1:nx, 1:ny] = -0.25 * (
        vorticity[0:nx - 1, 0:ny - 1] + vorticity[0:nx - 1, 1:ny] + vorticity[1:nx, 0:ny - 1] + vorticity[1:nx, 1:ny])
rhs[boundary == 1] = bdata

psi = LU(np.reshape(rhs.T, nxy_one))
psi = np.reshape(psi.T, (ny + 1, nx + 1)).T
psi[boundary == 1] = bdata

ux = 1 / dy * (psi[:, 1:ny + 1] - psi[:, 0:ny])
uy = -1 / dx * (psi[1:nx + 1, :] - psi[0:nx, :])
rk = np.array([[0, 1 / 2, 1 / 2, 1], [0, 1, 1, 1], [1, 2, 2, 1]]).T
ii = 0
uxs = np.zeros(ux.shape)  # uxs is our stored velocity from our previous time-steps
uys = np.zeros(uy.shape)  # used to compare our quasi steady state

#t = 0
resid = np.zeros(10000)  # Defining an array for storage of "residuals"
pos = 0
resid[pos] = np.mean((uxs - ux)) / np.mean(ux)

ux_max = np.max(ux)
uy_max = np.max(uy)
dt = np.min([dt_max_diffusion, 1 / (ux_max / dx + uy_max / dy), 10])

temp_value = 0
start_time_solver = time.time()
np.linspace(10, nx-10, 10).astype(int)
residual_position_x = np.linspace(10, nx-10, 10).astype(int) #np.array((int(nx/4),int(nx/2),int(3*nx/4)))
residual_position_y = np.linspace(ny-10, 10, 10).astype(int) #np.array((int(3*ny/4),int(ny/2),int(ny/4)))

ijk = 0
ux_storage = 1#np.zeros((nx+1,ny,10000))
uy_storage = 1#np.zeros((nx, ny+1, 10000))
vort_storage = np.zeros((nx, ny, 20001))
#ux_storage[:, :, 0] = ux
#uy_storage[:, :, 0] = uy
vort_storage[:,:,0] = vorticity
print('Time step is {} and steps {}'.format(dt, storage_steps))
minumum_timestep = 1*dt
t_storage = np.zeros(10000)
dt_min = 10
sub_storage = np.zeros(vort_storage.shape)

while t.current < RT + 10:
    ii += 1
    lt = time.time()
    irk = 0
    consump *= 0  # Reset consumption
    ux_max = np.max(ux)
    # print(ux_max)
    uy_max = np.max(uy)
    dt = np.min([dt_max_diffusion, 1 / (ux_max / dx + uy_max / dy), dt_min])
    # if ii % 30 == 0:
    #     uxs = ux
    #     uys = uy
    #     vort_storage[:, :, int(ijk)] = vorticity
    #     sub_storage[:, :, int(ijk)] = s.current
    #     t_storage[ijk] = t.current
    #     ijk += 1
    #     #print(ijk)


    while irk < 4:  # replaced with while loop to allow for adaptive timesteps
        if irk == 0:
            za.intermediate = za.current
            s.intermediate = s.current
            mox.intermediate = mox.current
            mred.intermediate = m_total - mox.intermediate
            vort1 = vorticity
            psi1 = psi
        else:
            za.update_intermediate(rk, irk, dt)
            s.update_intermediate(rk, irk, dt)
            mox.update_intermediate(rk, irk, dt)
            mred.intermediate = m_total - mox.intermediate
            vort1 = vorticity + rk[irk, 0] * dt * d_vorticity_dt1
            rhs[1:nx, 1:ny] = -0.25 * (
                    vort1[0:nx - 1, 0:ny - 1] + vort1[0:nx - 1, 1:ny] + vort1[1:nx, 0:ny - 1] + vort1[1:nx,
                                                                                                1:ny])
            rhs[boundary == 1] = bdata
            psi1 = LU(np.reshape(rhs.T, nxy_one))
            psi1 = np.reshape(psi1.T, (ny + 1, nx + 1)).T
            psi1[boundary == 1] = bdata

            # s.intermediate[0, round(int(1 / 2 * ny * 1 / (2 * baffle_pairs + 1) - in_out_points / 2)):round(
            #     int(1 / 2 * ny * 1 / (2 * baffle_pairs + 1) - in_out_points / 2 + 1 + in_out_points))] = s.influent

        if (mox.current + (dt / 6) * mox.ddt2 > m_total).any() or (
                mox.current + rk[irk, 0] * dt * mox.ddt1 > m_total).any():
            # If over estimating rk4 loop is reset with smaller timestep
            irk = 0
            dt *= 0.5
            continue

        ux = 1 / dy * (psi[:, 1:ny + 1] - psi[:, 0:ny])
        uy = -1 / dx * (psi[1:nx + 1, :] - psi[0:nx, :])
        fx[1:nx + 1, :] = ux[1:nx + 1, :] * vort1  # Initially assume everything is advected right
        fx[0, :] = ux[0, :] * vort1[0, :]
        reverse = (ux < 0)  # Check where advection is actually leftward
        reverse[nx, :] = 0
        fx[reverse] = ux[reverse] * vort1[reverse[0:nx, :]]  # Recalculate cells where advection is leftward
        fy[:, 1:ny + 1] = uy[:, 1:ny + 1] * vort1  # Initially assume everything is advected up
        fy[:, 0] = uy[:, 0] * vort1[:, 0]
        reverse = (uy < 0)  # Check where advection is actually downward
        reverse[:, ny] = 0
        fy[reverse] = uy[reverse] * vort1[reverse[:, 0:ny]]  # Recalculate cells where advection is downward
        trans_t = (1 / dx) * (fx[1:nx + 1, :] - fx[0:nx, :]) + (1 / dy) * (fy[:, 1:ny + 1] - fy[:, 0:ny])
        intermed = BB.dot(np.reshape(vort1.T, (nxy))) + CC.dot(np.reshape(psi1.T, nxy_one))
        intermed = np.reshape(intermed.T, (ny, nx)).T
        d_vorticity_dt1 = - trans_t + nu * intermed
        d_vorticity_dt2 = rk[irk, 1] * d_vorticity_dt2 + rk[irk, 2] * d_vorticity_dt1

        local_s = np.reshape(s.intermediate[bio_loc == 1], bio_number)
        substrate_mean_surface = anode_surface_sum(Ai * local_s, bio_loc) / A

        j, eta_act = current_density_inter(E, Rin, Rex, m_total, mred, mox, bio_loc, Ai, R, T, m, F, j0, full=BV_full)
        I_anode = anode_surface_sum_repeated(j * Ai, bio_loc)
        mediator_current_density_term = I_anode / A

        # local_s = np.reshape(s.intermediate[bio_loc == 1], (bio_number))
        za.update_growth_and_consumption(local_s, mox.intermediate)
        Xm.update_growth_and_consumption(local_s, mox.intermediate)
        diff_za = 0  # *Dx_bio.dot(vdata1[0,])
        diff_xm = 0  # *Dx_bio.dot(vdata1[1,]) # diffusion of X and Z
        za.first_timestep()
        za.second_timestep(rk, irk)
        mox.ddt1 = -Ym * za.consumption + pref * j / za.intermediate
        mox.second_timestep(rk, irk)

        s.calculate_advection(ux, uy, dx, dy)  # Advection is slowest process
        s.calculate_diffusion(Diffuse_s)  # diff is second slowest
        s.calculate_consumption(za, Xm, biofilm_location=bio_loc, convert_m2_l=convert_m2_l)  # rapid
        fluid_time = time.time()
        #s.consumption *= 0
        s.first_timestep()  # Timestepping is almost as slow as diff
        s.second_timestep(rk, irk)

        irk += 1  # move forward in rk4 loop
        if irk == 4 and (mox.current + (dt / 6) * mox.ddt2 > m_total).any():
            irk = 0
            dt *= 0.5
            # print('Loop restart')
            continue

    vorticity = vorticity + (dt / 6) * d_vorticity_dt2
    rhs[1:nx, 1:ny] = -0.25 * (
            vorticity[0:nx - 1, 0:ny - 1] + vorticity[0:nx - 1, 1:ny] + vorticity[1:nx, 0:ny - 1] + vorticity[1:nx,
                                                                                                    1:ny])
    rhs[boundary == 1] = bdata

    psi = LU(np.reshape(rhs.T, nxy_one))
    psi = np.reshape(psi.T, (ny + 1, nx + 1)).T
    psi[boundary == 1] = bdata

    # flux = flux_maximum / 4 * (np.cos(t / (1 * 12 * np.pi)) + 3)
    # flux = flux_maximum/2 * (np.cos(t / (30 * 12 * np.pi)) + 1)
    flux = flux_maximum / 2 * (np.cos(t.current / (30 * 12 * np.pi)) + 1)

    boundary_backup, psi_altered = boundary_init(14 * dy)
    psi_altered = system.influent_effluent_regions(baffle_pairs, baffle_length, 14 * dy, psi_altered,
                                                   boundary, flux)
    psi_altered = psi_altered[0]
    bdata_altered = psi_altered[boundary == 1]
    # psi[boundary == 1] = bdata_altered
    # rhs[boundary == 1] = bdata_altered

    ux = 1 / dy * (psi[:, 1:ny + 1] - psi[:, 0:ny])
    uy = -1 / dx * (psi[1:nx + 1, :] - psi[0:nx, :])
    bdata = 1 * bdata_altered

    za.update_current(dt)
    s.update_current(dt)
    mox.update_current(dt)
    mred.current = m_total - mox.current

    j, eta_act, eta_conc = current_density(E, Rin, Rex, m_total, mred, mox, bio_loc, Ai, R, T, m, F, j0, full=BV_full)

    I.current = np.sum(Ai * j)

    # surface_mean_biomass = anode_surface_sum_repeated(za.current * Ai, bio_loc) / A
    # Rin.current = Rin.minimum +(2000 - Rin.minimum)*np.exp(-0.0024*surface_mean_biomass)

    t.current += dt
    ux_max = ux.max()
    uy_max = uy.max()
    dt = min(dt_max, dt * 2, 1 / (ux_max / dx + uy_max / dy), (dx ** 2 * dy ** 2) / (
            2 * s.diffusion * (dx ** 2 + dy ** 2)))  # increase timestep up to 2 or double previous timestep
    # print('loop time ',t.current,time.time()-total_time)
    s.current = s.update_influent(baffle_pairs, in_out_points, ny)
    za.check_mass()
    mass_out = (ux[-1, (out_start - 1):out_start + in_out_points - 2] * (s.current[-1, (out_start - 1):out_start + in_out_points - 2] * Voli) / dx).sum()
    mass_in = (ux[0, (in_start - 1):in_start + in_out_points - 1] * s.current[0, (in_start - 1):in_start + in_out_points - 1] * Voli / dx).sum()
    #total_mass_in += mass_in*dt
    #total_mass_out += mass_out*dt
    #total_removed += -(s.consumption*Voli).sum()*dt
    #total_remaining = (s.current*Voli).sum()
    # if t.current > 4*TC_day and BV_full != True:
    #     BV_full = True

    # if (za.current > za.maximum).any():
    #     za.current[(za.current > za.maximum)] -= (za.current[(za.current > za.maximum)] - za.maximum)
    mass_flux = ((s.current[-2, out_start - 2:out_start + in_out_points] * ux[-1,out_start - 2:out_start + in_out_points]) * dy * Lz).sum()
    velocity_flux = ((ux[-1, out_start - 2:out_start + in_out_points]) * dy * Lz).sum()
    concentration_flux = mass_flux/velocity_flux



    if ii % storage_steps == 0 or ii == 1:  # round(t.now,2)%(k*60) == 0 : #Storage of data
        ij += 1
        za.storage[:, ij] = za.current
        Xm.storage[:, ij] = Xm.current
        mox.storage[:, ij] = mox.current
        mred.storage[:, ij] = mred.current
        I.storage[ij] = I.current
        t.storage[ij] = t.current
        s.storage[:, :, ij] = s.current
        s_blank[:, ij] = j#(-s.transported*Voli).sum() # (muz-Kda)*vdata[0,]
        increase = round((t.current - t.storage[ij - 1]) / (RT + 20) * 100, 1)
        pbar.update(round(increase, 1))
        s_blank1[ij] = concentration_flux

s.storage[s.storage>s.influent] = s.influent
# Plot_ani = True
# if Plot_ani == True:
#     fig, ax = plt.subplots(figsize=(16, 9))
#     cs = ax.contourf(xx, yy, s.storage[:, :, ij-1], cmap='coolwarm',levels = 50,vmin = 0,vmax = 100)
#     cs.cmap.set_over('#380e0b')
#     cs.changed()
#     fig.colorbar(cs)
#     import matplotlib.animation as animation
# ##### Home PC Directory ####
#     plt.rcParams['animation.ffmpeg_path'] = 'C:/Users/Jordan/Downloads/ffmpeg/ffmpeg/bin/ffmpeg'
# ##### Work Laptop Directory ###
#     #plt.rcParams['animation.ffmpeg_path'] = 'H:Downloads/ffmpeg/ffmpeg/bin/ffmpeg'
#
#
#
#     def animate(i):
#         ax.clear()
#         cs = ax.contourf(xx, yy, s.storage[:, :, i], cmap='coolwarm',levels = 50,vmin = 0,vmax = 100) #vmin = -0.1,vmax = 0.1
#         cs.cmap.set_over('#380e0b')
#         cs.changed()
#         ax.set_title('Time = {time:.2f} (Days)'.format(time = t.storage[i]/TC_day),fontsize = 20)
#         #fig.colorbar(cs)
#
#
#     interval = 1 / 30  # in seconds
#     ani = animation.FuncAnimation(fig, animate, ij, interval=interval * 1e-1, blit=False)
#     #plt.show()
#
# FFwriter=animation.FFMpegWriter(fps=60, extra_args=['-vcodec', 'libx264'])
# ani.save('1000mg_Pulse_timeseries_wide.mp4',writer=FFwriter)

def save_timeseries_classes(file_name, folder_names, substrate, biomass, mox, mred, I, t,s_blank,ij):
    """

    :param file_name: name for saved file
    :param folder_names: List of folder where file_name is to be saved
    :param substrate: substrate class
    :param biomass: biomass class
    :param mox: mox class
    :param mred: mred class
    :param j: current density array
    :param t: time class
    :return: will not return anything, instead saves file in location specified
    """

    file_path = Path(os.getcwd())
    for _ in np.arange(len(folder_names)):
        file_path = Path(file_path, folder_names[_])
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    file_path = Path(file_path, file_name)
    np.savez_compressed(file_path, Acetate=substrate.storage[...,0:ij]
                        , biomass=biomass.storage[...,0:ij]
                        , mox=mox.storage[...,0:ij]
        #                , mred=mred.storage[...,0:ij]
                        , current=I.storage[...,0:ij]
                        , time=t.storage[...,0:ij]
                        , current_density = s_blank[...,0:ij])

file_name = "100mg_Pulse_timeseries_low_maximum"
file_name = "100mg_Pulse_timeseries_low_maximum_5000_max_narrow"
save_timeseries_classes(file_name,[],s,za,mox,mred,I,t,s_blank,ij)



#
# save_data_classes(file_name,['Output','temp'],s,za,mox,mred,j,t)
# # save_data_classes('Single_substrate_cod_825',s,za,mox,mred,j,t)
# # save_data_classes('Single_substrate_cod_825',['Output','temp'],s,za,mox,mred,j,t)
#
# plot_positional_data(x, j, bio_loc, new_fig=True)
#
# print(current_density_inter(E, Rin, Rex, m_total, mred, mox, bio_loc, Ai, R, T, m, F, j0, full=False))
# current_density_inter(E, Rin, Rex, m_total, mred, mox, bio_loc, Ai, R, T, m, F, j0, full=True)
#
# print(dt, t.current)
# plot_positional_data(x, j, bio_loc, new_fig=True, side='Left',
#                      title='Current density using linear BV Eocv = {}'.format(E.current))
#
# plot_positional_data(x, eta_conc, bio_loc, side='left', new_fig=True,
#                      title='$\eta_\mathrm{conc} = M_\mathrm{total}/M_\mathrm{red}$')
#
# plt.figure(figsize=(14, 10))
# plt.subplot(221)
# plot_contour(xx, yy, s.current)
# plt.subplot(222)
# plot_positional_data(x, j, bio_loc, side='right', title='Positional current density A/m^2')
# plt.subplot(223)
# plot_time_series(t.storage[0:ij + 1] / TC_day, I.storage[0:ij + 1] / (20 * A), linelabel='Current density over time')
# plt.subplot(224)
# plot_positional_data(x, za.current, bio_loc, side='right', title='Positional biomass mg/m^2')
# # save_figure('E_applied_0_9_Rin_12_Rext_1_hrt_48')
#
#
# # (ux[0,(in_start-1):in_start+in_out_points-1]*(s.influent*in_out_points*Voli)/dx*dt).sum()
# # (ux[-1,(out_start-1):out_start+in_out_points-2]*(s.current[-1,(out_start-1):out_start+in_out_points-2]*Voli)/dx*dt).sum()
#
#
# # (za.consumption*za.current).sum()
# # (ux[0,(in_start-1):in_start+in_out_points-1]*(s.influent*in_out_points*Voli)/dx).sum()
# # (ux[-1,(out_start-1):out_start+in_out_points-2]*(s.current[-1,(out_start-1):out_start+in_out_points-2]*Voli)/dx).sum()
# # (ux[0,(in_start-1):in_start+in_out_points-1]*s.current[0,(in_start-1):in_start+in_out_points-1]*Voli).sum()
# # (ux[-1,(out_start-1):out_start+in_out_points-2]*(s.current[-1,(out_start-1):out_start+in_out_points-2]*Voli)).sum()
# mass_out = (ux[-1,(out_start-1):out_start+in_out_points-2]*(s.current[-1,(out_start-1):out_start+in_out_points-2]*Voli)/dx).sum()
# mass_in = (ux[0,(in_start-1):in_start+in_out_points-1]*s.current[0,(in_start-1):in_start+in_out_points-1]*Voli/dx).sum()
#
#
# Plot_ani = True
# if Plot_ani == True:
#     fig, ax = plt.subplots(figsize=(16, 9))
#     cs = ax.contourf(xx, yy, vorticity[:, :, 0], cmap='coolwarm',levels = 50)#,vmin = -0.15,vmax = 0.15)
#     cs.cmap.set_over('#380e0b')
#     cs.changed()
#     fig.colorbar(cs)
#     import matplotlib.animation as animation
# ##### Home PC Directory ####
#     plt.rcParams['animation.ffmpeg_path'] = 'C:/Users/Jordan/Downloads/ffmpeg/ffmpeg/bin/ffmpeg'
# ##### Work Laptop Directory ###
#     #plt.rcParams['animation.ffmpeg_path'] = 'H:Downloads/ffmpeg/ffmpeg/bin/ffmpeg'
#
#
#
#     def animate(i):
#         ax.clear()
#         cs = ax.contourf(xx, yy, vorticity[:, :, i], cmap='coolwarm',levels = 50,vmin = -0.15,vmax = 0.15) #vmin = -0.1,vmax = 0.1
#         cs.cmap.set_over('#380e0b')
#         cs.changed()
#         ax.set_title('Time = {}'.format(t[i]),fontsize = 20)
#         #fig.colorbar(cs)
#
#
#     interval = 1 / 30  # in seconds
#     ani = animation.FuncAnimation(fig, animate, 1000, interval=interval * 1e-1, blit=False)
#     #plt.show()
#
# FFwriter=animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
# ani.save('Vorticity_Fluctuating_6hr_0_period_2hrs.mp4',writer=FFwriter)
#
#
# def save_timeseries_classes(file_name, folder_names, substrate, biomass, mox, mred, I, t,s_blank,ij):
#     """
#
#     :param file_name: name for saved file
#     :param folder_names: List of folder where file_name is to be saved
#     :param substrate: substrate class
#     :param biomass: biomass class
#     :param mox: mox class
#     :param mred: mred class
#     :param j: current density array
#     :param t: time class
#     :return: will not return anything, instead saves file in location specified
#     """
#
#     file_path = Path(os.getcwd())
#     for _ in np.arange(len(folder_names)):
#         file_path = Path(file_path, folder_names[_])
#     if not os.path.exists(file_path):
#         os.makedirs(file_path)
#     file_path = Path(file_path, file_name)
#     np.savez_compressed(file_path, Acetate=substrate.storage[...,0:ij]
#                         , biomass=biomass.storage[...,0:ij]
#                         , mox=mox.storage[...,0:ij]
#         #                , mred=mred.storage[...,0:ij]
#                         , current=I.storage[...,0:ij]
#                         , time=t.storage[...,0:ij]
#                         , current_density = s_blank[...,0:ij])
#
# file_name = "100mg_Pulse_timeseries"/
# save_timeseries_classes(file_name,[],s,za,mox,mred,I,t)