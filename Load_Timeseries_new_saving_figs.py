# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 10:02:47 2022

@author: Jordan
"""


import numpy as np 
from Fluid_classes import *
from Plotting_functions import *
import os
from pathlib import Path

def loading_data(file_name, Folder_names ='Empty'):
    """
    :param file_name: Name of loading data file
    :param Folder_names: List of folders where file_name is located
    :return: loaded substrate, biomass, mox, mred and current density
    """

    if Folder_names == 'Empty':
        loaded_data = np.load(file_name)
        substrate = loaded_data['Acetate']
        biomass = loaded_data['biomass']
        mox = loaded_data['mox']
        #mred = loaded_data['mred']
        j = loaded_data['current_density']
        I = loaded_data['current']
        return substrate,biomass,mox,mred,j
    else:
        file_path = Path(os.getcwd())
        for _ in np.arange(len(Folder_names)):
            file_path = Path(file_path, Folder_names[_])
        file_path = Path(file_path, file_name)
        loaded_data = np.load(file_path)
        substrate = loaded_data['Acetate']
        biomass = loaded_data['biomass']
        mox = loaded_data['mox']
        #mred = loaded_data['mred']
        j = loaded_data['current_density']
        I = loaded_data['current']
        t = loaded_data['time']
        return substrate, biomass, mox, I, j,t
    
    
def plot_streamlines(xmesh,ymesh,ux_data,uy_data,title = '',xlab = 'x',ylab= 'y',new_fig = False,stream_density = 3,cmapping = 'Greys'):
    if new_fig == True:
        plt.figure(figsize = (11,14))#figsize = (10,10*ymesh.max()/xmesh.max())) #figsize = (10,12)
    
    ux_i = 0.5*(ux_data[0:ux_data.shape[0]-1,:] + ux_data[1:ux_data.shape[0],:])
    uy_i = 0.5*(uy_data[:,0:uy_data.shape[1]-1] + uy_data[:,1:uy_data.shape[1]])
    speed = np.sqrt((ux_i.T)**2 +(uy_i.T)**2 )
    lw = 10*speed / speed.max()            
    plt.streamplot(xmesh.T, ymesh.T, ux_i.T, uy_i.T,color = speed, density=stream_density,linewidth=1.5,arrowstyle = '-|>')   #,cmap=cmapping
    cb = plt.colorbar(fraction=0.062, pad=0.04)
    cb.ax.locator_params(nbins=6)
    cb.ax.tick_params(labelsize=28)
    cb.ax.set_yticklabels(['$2.0 \cdot 10^{-6}$','$4.0 \cdot 10^{-6}$','$6.0 \cdot 10^{-6}$','$8.0 \cdot 10^{-6}$','$10 \cdot 10^{-6}$'])#,'$1.2 \cdot 10^{-3}$'])
    cb.ax.set_ylabel('m/s', rotation=270,fontsize = 28,y = 0.5)
    set_axis_titles(title, xlab, ylab, xmesh.max(), ymesh.max(),y_rotation = 0)
    plt.axis('scaled')
    plt.xlim([0,xmesh.max()])
    plt.ylim([0,ymesh.max()])
    plt.xticks(np.arange(0.1, xmesh.max(), 0.1))
    plt.yticks(np.arange(0, ymesh.max(), 0.1))
    plt.tight_layout()#h_pad=1)
    


effluent_substrate = np.zeros(6)
average_cd = np.zeros(6)
loaded_hrt = np.zeros(6)
loaded_baffle_pairs = np.zeros(6)

nx = 300
ny = 300
nz = 1
# Size
Lx = 1 # length in m
Ly = 1 # length in m
Lz = 1 # length in m
# HRT
hrt = 6  # hrt set in any given time units TC is the conversion to seconds
hrt *= TC_hour  # The hydraulic retention time converted into seconds
#Baffle param
baffle_length = 91/100 ##### This is the fraction of the tank a baffle takes up in x
baffle_pairs = 2 ## 2 for wide, 6 for narrow  #### Number of baffle pairs (RHS+LHS) = 1 pair.
# Runtime for reactor system
RT  = 7
RT *= TC_day
dt_max = 4
k = 20 # Store data every k min, 1 is every min, 2 every 2 min
D = 1/hrt
stream = False
anode_side = 'all' # can either be 'all', 'influent', 'effluent'
Vol = int(Lx*Ly*Lz*1e3)


dx = Lx/nx
dy = Ly/ny
x = np.linspace(0,Lx,nx).T
y = np.linspace(0,Ly,ny).T
[yy,xx] = np.meshgrid(np.linspace(dy/2,Ly-dy/2,ny),np.linspace(dx/2,Lx-dx/2,nx))

system = domain(Lx, Ly, Lz,nx = nx,ny = ny)


flux = (Lx*Ly)/hrt # The "area flux" through the system
nxy = nx*ny
nxy_one = (nx+1)*(ny+1)
psi = np.zeros((nx+1,ny+1))  # This is the stream function
boundary = np.zeros((nx+1,ny+1))  # We set this to 1 on all boundary points
boundary[0,:] = 1
boundary[-1,:] = 1
boundary[:,0] = 1
boundary[:,-1] = 1

edges = boundary
psi[0,0:ny+3] = flux
psi[:,-1] = flux
psi,boundary,in_out_points,in_start,out_start = system.influent_effluent_regions(baffle_pairs, baffle_length, dy*14, psi, boundary, flux)
bdata = psi[boundary == 1]
external = np.zeros(boundary.shape)
external[0,:] = 1
external[-1,:] = 1
external[:,0] = 1
external[:,-1] = 1
internal = boundary - external
bio_loc = np.zeros(boundary.shape)
bio_loc[:,0:ny] += internal[:,1:ny+1]
bio_loc[:,1:ny+1] += internal[:,0:ny]
bio_loc = bio_loc[0:nx,0:ny]
if baffle_length==0 and baffle_pairs == 0:
    bio_loc[1:-1,-2] = 1
    bio_loc[1:-1,1] = 1
file_a = 'hrt' + str(hrt).replace('.','_') + '_nx' +str(nx) + '_ny' +str(ny)+ '_Lx' +str(Lx) + '_Ly' +str(Ly)+'_pairs'+str(baffle_pairs)+'_width'+str(np.round(baffle_length,decimals=1)).replace('.','_')+'.csv'
file_x = 'Ux_' + file_a
file_y = 'Uy_' + file_a
data_folder = Path(os.getcwd(),"Output","Velocity")
#ux = np.genfromtxt(data_folder / file_x, delimiter=',')
#uy = np.genfromtxt(data_folder / file_y, delimiter=',')
anode_side = 'all'
pos = baffle_pairs - 1

#data_file_name =  'Y_22_s_500_RT_7_days_HRT_48.npz'#.format(int(hrt/TC_hour)) # 'Y_22_s_500_RT_7_days_HRT_12_baffles_6.npz'
#data_file_name =  'Low_Substrate_4_baffles_4Days.npz'#.format(int(hrt/TC_hour))
#WIDE DATASET 
data_file_name =  '100mg_Pulse_timeseries_wide.npz'#.format(int(hrt/TC_hour))
#NARROW DATASET 
#data_file_name =  '100mg_Pulse_timeseries_narrow.npz'#.format(int(hrt/TC_hour))500mg_Pulse_timeseries_inter
#Intermediate DATASET 
#data_file_name =  '500mg_Pulse_timeseries_inter.npz'
#s,za,mox,I,j,t = loading_data(data_file_name,[])

Ai = dx * Lz
A = baffle_length * nx * Ai


flux_maximum = 1*flux

#%%
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

#flux_maximum = 1*flux
mass_flux = np.zeros((t.shape[0]))
velocity_flux = np.zeros((t.shape[0]))
stored_flux = np.zeros((t.shape[0]))

for ijkl in np.arange(t.shape[0]):
    flux = flux_maximum / 2 * (np.cos(t[ijkl] / (30 * 12 * np.pi)) + 1)
    

    boundary_backup, psi_altered = boundary_init(14 * dy)
    psi_altered = system.influent_effluent_regions(baffle_pairs, baffle_length, 14 * dy, psi_altered,
                                               boundary, flux)
    psi_altered = psi_altered[0]
    bdata_altered = psi_altered[boundary == 1]
# psi[boundary == 1] = bdata_altered
# rhs[boundary == 1] = bdata_altered

    ux = 1 / dy * (psi_altered[:, 1:ny + 1] - psi_altered[:, 0:ny])
    bdata = 1 * bdata_altered
    mass_flux[ijkl] = ((s[-2, out_start-2:out_start + in_out_points,ijkl]*ux[-1, out_start-2:out_start + in_out_points])*dy*Lz).sum()
    velocity_flux[ijkl] = ((ux[-1, out_start-2:out_start + in_out_points])*dy*Lz).sum()
    stored_flux[ijkl] = flux

plt.plot(t/(60**2*24),(s[...,1].max() - mass_flux/velocity_flux)/s[...,1].max()*100,label='1000mg Wide')

#%%
data_file_name =  '4baffles_100influent_pulse.npz' #'100mg_Pulse_timeseries_low_maximum_5000_max.npz'
#100mg_Pulse_timeseries_low_maximum_5000_max_wide
s,za,mox,I,j,t = loading_data(data_file_name,[])

s[s>100] = 100

baffle_pairs = 6
#print(hello_world)
f = plt.figure(1,figsize(24,12))
f, (ax3, ax4,ax5,ax6) = plt.subplots(1, 4)
for i in np.arange(0,6375,2):

    ax3.clear()
    ax3.set_ylim(-0.1,1.5)
    ax3.set_xlim(-0.1,4.1)
    ax3.plot(t[0:i]/60/60/24,I[0:i]/(baffle_pairs*4*A),linewidth = 3)
    ax3.set_ylabel('Current Density (A/$m^2$)',fontsize = 18)
    ax3.set_xlabel('Time (Days)',fontsize = 18)
    ax3.tick_params(axis='both', labelsize=16)
    #ax3.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    ax3.set_title('Average Current Density',fontsize = 20)
    
    ax4.clear()
    ## Input function for plotting positional data
    biofilm_side = int(x.size / 2)
    positional_layout = np.zeros(bio_loc.shape)
    positional_layout[bio_loc == 1] = j[...,i]
    y_positions = np.arange(bio_loc.shape[1])[bio_loc[biofilm_side, :] == 1]
    for position in np.arange(y_positions.size):
        plotting_data = positional_layout[:, y_positions[position]]
        plotting_data[plotting_data == 0] = np.nan  # Marking zero data as zero to clean up plotting
        ax4.plot(x, plotting_data, '-', label=position + 1, linewidth=3)  # position+1 linelabel
    ax4.set_ylim(-0.1,3)
    ax4.set_xlim(-0.1,1.1)
    #ax4.plot(j[...,i],'.')
    ax4.set_ylabel('Current Density (A/$m^2$)',fontsize = 18)
    ax4.set_xlabel('Position',fontsize = 18)
    #ax4.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    ax4.set_title('Positional Current Density',fontsize = 20)
    ax4.tick_params(axis='both', labelsize=16)
    
    ax5.clear()
    ## Input function for plotting positional data
    ax5.set_ylim(-100,5100)
    ax5.set_xlim(-0.1,1.1)
    #ax5.plot(za[...,i],'.g')
    biofilm_side = int(x.size / 2)
    positional_layout = np.zeros(bio_loc.shape)
    positional_layout[bio_loc == 1] = za[...,i]
    y_positions = np.arange(bio_loc.shape[1])[bio_loc[biofilm_side, :] == 1]
    for position in np.arange(y_positions.size):
        plotting_data = positional_layout[:, y_positions[position]]
        plotting_data[plotting_data == 0] = np.nan  # Marking zero data as zero to clean up plotting
        ax5.plot(x, plotting_data, '-', label=position + 1, linewidth=3)  # position+1 linelabel
    ax5.set_ylabel('Biomass (mg/$m^2$)',fontsize = 18)
    ax5.set_xlabel('Position',fontsize = 18)
    #ax5.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    ax5.set_title('Positional Biomass',fontsize = 20)
    ax5.tick_params(axis='both', labelsize=16)
    
    # ax3.clear()
    # ax4.clear()
    # ax5.clear()
    ax6.clear()
    cs = ax6.contourf(xx, yy, s[:, :, i], cmap='coolwarm',levels = 50,vmin = 0,vmax = 100) #vmin = -0.1,vmax = 0.1
    ax6.set_ylabel('y',fontsize = 18,rotation = 0)
    ax6.set_xlabel('x',fontsize = 18)
    cs.cmap.set_over('#380e0b')
    cs.changed()
    ax6.set_title('Substrate Contour at t = {time:.2f} Days'.format(time = t[i]/TC_day),fontsize = 20)
    plt.tight_layout()
    ax6.tick_params(axis='both', labelsize=16)
    save_figure('hrt_vs_cd_substrate_removal_Y_100_'+str(i))
    
#%%
s[s>100] = 100

baffle_pairs = 2
f = plt.figure(1,figsize(10,8))
f, ax = plt.subplots(1, 1)

cs = ax.contourf(xx, yy, s[:, :, 100], cmap='coolwarm',levels = 50,vmin = 0,vmax = 100) #vmin = -0.1,vmax = 0.1
cbar = plt.colorbar(cs)
cbar.set_label('Concentration (mg/L)',fontsize = 36)
cbar.ax.tick_params(labelsize = 30)
ax.clear()


#for i in np.arange(0,za.shape[1],100):
for i in [4642,4652,4664,4678,4693]:

    
    # ax3.clear()
    # ax4.clear()
    # ax5.clear()
    ax.clear()
    cs = ax.contourf(xx, yy, s[:, :, i], cmap='coolwarm',levels = 50,vmin = 0,vmax = 100) #vmin = -0.1,vmax = 0.1
    ax.set_ylabel('y',fontsize = 36,rotation = 0)
    ax.set_xlabel('x',fontsize = 36)
    cs.cmap.set_over('#380e0b')
    cs.changed()
    ax.set_title('t = {time:.2f} Hours'.format(time = t[i]/(60**2)),fontsize = 36)
    plt.tight_layout()
    ax.tick_params(axis='both', labelsize=30)
    plt.tight_layout()
    save_figure('substrate_contour_Y_100_'+str(i))
    


#%%
baffle_pairs = 2
f = plt.figure(1,figsize(12,10))
f, ax3 = plt.subplots(1, 1)
ax3.clear()
ax3.set_ylim(-0.1,1.5)
ax3.set_xlim(-0.1,4.1)
ax3.plot(t[:]/60/60/24,I[:]/(baffle_pairs*4*A),'b',linewidth = 4,)
ax3.set_ylabel('Current Density (A/$m^2$)',fontsize = 36)
ax3.set_xlabel('Time (Days)',fontsize = 36)
ax3.tick_params(axis='both', labelsize=16)
#ax3.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
#ax3.set_title('Average Current Density',fontsize = 36)
import pandas as pd
I_df = pd.DataFrame(I/(baffle_pairs*4*A))
ax3.plot(t[:]/60/60/24,I_df.rolling(100).mean(),'-.k',linewidth = 6)

cont_dataset = np.genfromtxt('Timeseries_data_wide_100_cont.txt',delimiter=",")
cont_I = np.trim_zeros(cont_dataset[0])
cont_t = np.trim_zeros(cont_dataset[1])
plt.plot((cont_t/60/60/24)[cont_t/60/60/24 < 4],(cont_I/(baffle_pairs*4*A))[cont_t/60/60/24 < 4],'--r',linewidth = 6)
ax3.tick_params(axis='both', labelsize=30)
plt.tight_layout()
ax3.legend(['Pulsing','Pulsing (Moving Average)','Continuous'],fontsize = 30)

#%%
def mean_effluent(s, in_out_points, out_start):
    temp = np.mean(s[-2, out_start:out_start + in_out_points])
    return temp

import pandas as pd



effluent = np.zeros(s.shape[-1])

for ijk in np.arange(s.shape[-1]):
    effluent[ijk] = mean_effluent(s[:,:,ijk],in_out_points,out_start)
    #effluent[ijk] = np.sum(s[...,ijk]*dx*dy*Lz)


#plt.plot(t[:]/60/60/24,(100 - effluent)/100*100)
#plt.plot((cont_t/60/60/24)[cont_t/60/60/24 < 4],((100 - cont_effluent)/100*100)[cont_t/60/60/24 < 4])

f = plt.figure(1,figsize(12,10))
f, ax3 = plt.subplots(1, 1)
ax3.clear()
#ax3.set_ylim(-0.1,1.5)
#ax3.set_xlim(-0.1,4.1)
ax3.plot(t[:]/60/60/24,(100 - effluent)/100*100,'b',linewidth = 4,)
ax3.set_ylabel('Substrate Removal (%)',fontsize = 36)
ax3.set_xlabel('Time (Days)',fontsize = 36)
#ax3.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
#ax3.set_title('Average Current Density',fontsize = 36)
effluent_df = pd.DataFrame(effluent)
ax3.plot(t[:]/60/60/24,(100 - effluent_df.rolling(100).mean())/1,'-.k',linewidth = 6)

cont_effluent = np.genfromtxt('effluent_100_wide_cont.txt',delimiter=",")
plt.plot((cont_t/60/60/24)[cont_t/60/60/24 < 4],((100 - cont_effluent)/100*100)[cont_t/60/60/24 < 4],'--r',linewidth = 6)
ax3.tick_params(axis='both', labelsize=30)
plt.tight_layout()
ax3.legend(['Pulsing','Pulsing (Moving Average)','Continuous'],fontsize = 30)
#%%

f = plt.figure(1,figsize(12,10))
f, ax3 = plt.subplots(1, 1)
ax3.clear()
#ax3.set_ylim(-0.1,1.5)
#ax3.set_xlim(-0.1,4.1)
ax3.plot(t[:]/60/60/24,(100 - effluent)/100*100,'b',linewidth = 4,)
ax3.set_ylabel('Substrate Removal (%)',fontsize = 36)
ax3.set_xlabel('Time (Days)',fontsize = 36)
#ax3.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
#ax3.set_title('Average Current Density',fontsize = 36)
effluent_df = pd.DataFrame(effluent)
#ax3.plot(t[:]/60/60/24,(100 - effluent_df.rolling(100).mean())/100*100,'-.k',linewidth = 6)
ax3.tick_params(axis='both', labelsize=30)
plt.tight_layout()
ax3.set_title('Volume Integral Over Full Domain',fontsize = 36)
plt.tight_layout()
ax3.legend(['Pulsing','Pulsing (Moving Average)','Continuous'],fontsize = 30)



#%%

cont_substrate = np.genfromtxt('substrate_100_wide_cont.txt',delimiter=",")

baffle_pairs = 2
f = plt.figure(1,figsize(10,8))
f, ax = plt.subplots(1, 1)

cs = ax.contourf(xx, yy, cont_substrate, cmap='coolwarm',levels = 50,vmin = 0,vmax = 100) #vmin = -0.1,vmax = 0.1
cbar = plt.colorbar(cs)
cbar.set_label('Concentration (mg/L)',fontsize = 36)
cbar.ax.tick_params(labelsize = 30)
ax.clear()

ax.clear()
cs = ax.contourf(xx, yy, cont_substrate, cmap='coolwarm',levels = 50,vmin = 0,vmax = 100) #vmin = -0.1,vmax = 0.1
ax.set_ylabel('y',fontsize = 36,rotation = 0)
ax.set_xlabel('x',fontsize = 36)
cs.cmap.set_over('#380e0b')
cs.changed()
ax.set_title('12 Hour HRT',fontsize = 36)
plt.tight_layout()
ax.tick_params(axis='both', labelsize=30)
plt.tight_layout()
# save_figure('substrate_contour_12hr_HRT')





#%%
data_file_name =  '4baffles_100influent_pulse.npz' #'100mg_Pulse_timeseries_low_maximum_5000_max.npz'
#100mg_Pulse_timeseries_low_maximum_5000_max_wide
s,za,mox,I,j,t = loading_data(data_file_name,[])
print(s.max(),za.max())

data_file_cont = '4baffles_100influent_cont.npz'
s1,za1,mox1,I1,j1,t1 = loading_data(data_file_cont,[])

dt = t[1:] - t[0:-1] # 1 Shorter than dataset
un_inter_j = I/(4*baffle_pairs*A)
inter_j_time = np.zeros(un_inter_j.shape)

store_gap = 200

for ikm in np.arange(store_gap+1,I.shape[0]):
    temp_sum = un_inter_j[ikm-store_gap:ikm]
    temp_dt = dt[ikm-store_gap-1:ikm -1]
    temp_intergral = np.sum(temp_sum *temp_dt) / (t[ikm] - t[ikm-store_gap])
    inter_j_time[ikm] = temp_intergral
    #del temp_sum, temp_dt, temp_intergral
    

f = plt.figure(1,figsize(12,10))
f, ax3 = plt.subplots(1, 1)  
    
    
ax3.plot(t[:]/60/60/24,un_inter_j,'b',linewidth = 4)
ax3.plot(t[:][store_gap+10:]/60/60/24,inter_j_time[store_gap+10:],'-.k',linewidth = 6)
ax3.plot((t1/60/60/24)[t1/60/60/24 < 4],(I1/(baffle_pairs*4*A))[t1/60/60/24 < 4],'--r',linewidth = 6)
ax3.set_ylabel('Current Density (A/$m^2$)',fontsize = 36)
ax3.set_xlabel('Time (Days)',fontsize = 36)
ax3.tick_params(axis='both', labelsize=30)
ax3.legend(['Pulsing','Pulsing (Integral Average)','Continuous'],fontsize = 30)
plt.tight_layout()


flux_maximum = 1*flux
effluent_updated = np.zeros(t.shape)
effluent_updated1 = np.zeros(t1.shape)

mass_flux_store = np.zeros(t.shape)
velocity_flux_store = np.zeros(t.shape)

for i in np.arange(t.shape[0]):
    flux = flux_maximum / 2 * (np.cos(t[i] / (30 * 12 * np.pi)) + 1)
    boundary_backup, psi_altered = boundary_init(14 * dy)
    psi_altered = system.influent_effluent_regions(baffle_pairs, baffle_length, 14 * dy, psi_altered,
                                                   boundary, flux)
    psi_altered = psi_altered[0]
    bdata_altered = psi_altered[boundary == 1]
    # psi[boundary == 1] = bdata_altered
    # rhs[boundary == 1] = bdata_altered

    ux = 1 / dy * (psi[:, 1:ny + 1] - psi[:, 0:ny])

    mass_flux = ((s[...,i][-2, out_start - 2:out_start + in_out_points] * ux[-1,out_start - 2:out_start + in_out_points]) * dy * Lz).sum()
    velocity_flux = ((ux[-1, out_start - 2:out_start + in_out_points]) * dy * Lz).sum()
    concentration_flux = mass_flux/velocity_flux
    mass_flux_store[i] = mass_flux
    velocity_flux_store[i] = velocity_flux
    effluent_updated[i] = concentration_flux





dt = t[1:] - t[0:-1] # 1 Shorter than dataset
un_inter_s = effluent_updated
inter_s_time = np.zeros(effluent_updated.shape)
inter_s_time2 = np.zeros(effluent_updated.shape)

store_gap = 200

for ikm in np.arange(store_gap+1,I.shape[0]):
    temp_sum = un_inter_s[ikm-store_gap:ikm]
    
    temp_sum_con =  mass_flux_store[ikm-store_gap:ikm]
    temp_sum_vel = velocity_flux_store[ikm-store_gap:ikm]
    
    
    temp_dt = dt[ikm-store_gap-1:ikm -1]
    temp_intergral = np.sum(temp_sum *temp_dt) / (t[ikm] - t[ikm-store_gap])
    inter_s_time[ikm] = temp_intergral
    inter_s_time2[ikm] = np.sum(temp_sum_con*temp_dt)/np.sum(temp_sum_vel*temp_dt)
    #del temp_sum, temp_dt, temp_intergral
    

flux = flux_maximum / 2 
boundary_backup, psi_altered = boundary_init(14 * dy)
psi_altered = system.influent_effluent_regions(baffle_pairs, baffle_length, 14 * dy, psi_altered,
                                                   boundary, flux)
psi_altered = psi_altered[0]
bdata_altered = psi_altered[boundary == 1]
ux = 1 / dy * (psi[:, 1:ny + 1] - psi[:, 0:ny])
velocity_flux = ((ux[-1, out_start - 2:out_start + in_out_points]) * dy * Lz).sum()
for _ in np.arange(s1.shape[2]):
    
    mass_flux = ((s1[...,_][-2, out_start - 2:out_start + in_out_points] * ux[-1,out_start - 2:out_start + in_out_points]) * dy * Lz).sum()
    
    concentration_flux = mass_flux/velocity_flux
    effluent_updated1[_] = concentration_flux


temp_sum = ((inter_s_time2))[1200:-2]
temp_dt = dt[1200:-1]
np.sum(temp_sum *temp_dt) / (t[-2] - t[1200])




f = plt.figure(1,figsize(12,10))
f, ax3 = plt.subplots(1, 1)  
    
    
ax3.plot(t[:]/60/60/24,(s1.max()- un_inter_s)/s1.max()*100,'b',linewidth = 4)
ax3.plot(t[:][store_gap+10:]/60/60/24,(s1.max()-inter_s_time2[store_gap+10:])/s1.max()*100,'-.k',linewidth = 6)
ax3.legend(['1','2'])
ax3.plot((t1/60/60/24)[t1/60/60/24 < 4],((s1.max()-effluent_updated1)/s1.max()*100)[t1/60/60/24 < 4],'--r',linewidth = 6)
ax3.set_ylabel('Substrate Removal (%)',fontsize = 36)
ax3.set_xlabel('Time (Days)',fontsize = 36)
ax3.tick_params(axis='both', labelsize=30)
ax3.legend(['Pulsing','Pulsing (Integral Average)','Continuous'],fontsize = 30)
plt.tight_layout()


