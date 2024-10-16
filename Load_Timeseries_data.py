# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 11:16:43 2022

@author: Jordan
"""



# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 21:38:20 2020

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
baffle_pairs = 6 #### Number of baffle pairs (RHS+LHS) = 1 pair.
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
# data_file_name =  'Low_Substrate_4_baffles_4Days.npz'#.format(int(hrt/TC_hour))
#WIDE DATASET 
data_file_name =  '1000mg_Pulse_timeseries_wide.npz'#.format(int(hrt/TC_hour))
#NARROW DATASET 
# data_file_name =  '100mg_Pulse_timeseries_narrow.npz'#.format(int(hrt/TC_hour))
s,za,mox,I,j,t = loading_data(data_file_name,[])

Ai = dx * Lz
A = baffle_length * nx * Ai

#plot_contour(xx,yy,s,new_fig=True,ux_data = ux,uy_data = uy,scatter_data = za,scatter_positions= bio_loc,title = '12 Cartridges',stream_density = 4) #'Hrt 2 Hour' '12 Cartridges' 'Hrt 1 Hour'
#save_figure('hrt_48_substrate_contour_biomass')
#save_figure('hrt_12_cartridges_12')

# plot_contour(xx,yy,s,new_fig=True,scatter_data = za,scatter_positions= bio_loc,title = '',stream_density = 4) #'Hrt 2 Hour' '12 Cartridges' 'Hrt 1 Hour'
#plot_streamlines(xx,yy,ux,uy,new_fig=True,title = 'HRT 240 Hours',stream_density = 5)

#%%
import matplotlib.animation as animation
#from matplotlib.animation import FuncAnimation 
fig, ax = plt.subplots()
x_data=[]
y_data=[]

ax.set_ylim(-0.1,1.5)
ax.set_xlim(-0.1,4.1)
line, = ax.plot(0,0)

def animation_frame(i):
    global t, I,A
    x_new = t[int(i)]/60/60/24
    y_new = I[int(i)]/(8*A)
    x_data.append(x_new)
    y_data.append(y_new)
    
    line.set_xdata(x_data)
    line.set_ydata(y_data)
    ax.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    return line,

ani = animation.FuncAnimation(fig,func=animation_frame,frames = len(t),interval = 0.05,repeat=False)
plt.show()
plt.rcParams['animation.ffmpeg_path'] = 'C:/Users/Jordan/Downloads/ffmpeg/ffmpeg/bin/ffmpeg'
FFwriter=animation.FFMpegWriter(fps=60, extra_args=['-vcodec', 'libx264'])
#ani.save('Current_Density_Animation_100mg_wide.mp4',writer=FFwriter)

#%%
import matplotlib.animation as animation
fig1, ax1 = plt.subplots()
ax1.set_ylim(-0.1,3)
ax1.set_xlim(-0.1,2180)

def ani_frame(i):
    global j
    ax1.clear()
    ax1.set_ylim(-0.1,3)
    ax1.set_xlim(-0.1,2180)
    ax1.plot(j[...,i],'.')
    ax1.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    
ani = animation.FuncAnimation(fig1,func=ani_frame,frames = j.shape[1],interval = 2,repeat=True)
plt.show()
ani.save('Current_Density_Positional_Dots_Animation_1000mg_wide.mp4',writer=FFwriter)


#%%

f, (ax3, ax4) = plt.subplots(1, 2)
def ani_frame1(i):
    global t, I,A,j
    ax3.clear()
    ax3.set_ylim(-0.1,1.3)
    ax3.set_xlim(-0.1,4.1)
    ax3.plot(t[0:i]/60/60/24,I[0:i]/(8*A))
    ax3.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    
    ax4.clear()
    ax4.set_ylim(-0.1,3)
    ax4.set_xlim(-0.1,2180)
    ax4.plot(j[...,i],'.')
    ax4.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    

ani = animation.FuncAnimation(f,func=ani_frame1,frames = j.shape[1],interval = 0.5,repeat=False)
#ani.save('Current_Density_combined_Animation_100mg_wide.mp4',writer=FFwriter)

#%%

fig1, ax1 = plt.subplots()
ax1.set_ylim(-0.1,3)
ax1.set_xlim(-0.1,2180)

def ani_frame(i):
    global za
    ax1.clear()
    ax1.set_ylim(4000,7100)
    ax1.set_xlim(-0.1,2180)
    ax1.plot(za[...,i],'.g')
    ax1.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    
ani = animation.FuncAnimation(fig1,func=ani_frame,frames = za.shape[1],interval = 2,repeat=True)
plt.show()
ani.save('Biomass_Positional_Dots_Animation_1000mg_wide.mp4',writer=FFwriter)

#%%
f = plt.figure(figsize(24,12))
f, (ax3, ax4,ax5,ax6) = plt.subplots(1, 4)
def ani_frame1(i):
    global t, I,A,j,za,s,xx,yy
    ax3.clear()
    ax3.set_ylim(-0.1,1.5)
    ax3.set_xlim(-0.1,4.1)
    ax3.plot(t[0:i]/60/60/24,I[0:i]/(8*A))
    #ax3.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    ax3.set_title('Average Current Density',fontsize = 20)
    
    ax4.clear()
    ax4.set_ylim(-0.1,3)
    ax4.set_xlim(-0.1,2180)
    ax4.plot(j[...,i],'.')
    #ax4.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    ax4.set_title('Positional Current Density',fontsize = 20)
    
    ax5.clear()
    ax5.set_ylim(4000,7100)
    ax5.set_xlim(-0.1,2180)
    ax5.plot(za[...,i],'.g')
    #ax5.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    ax5.set_title('Positional Biomass',fontsize = 20)
    
    ax.clear()
    cs = ax6.contourf(xx, yy, s[:, :, i], cmap='coolwarm',levels = 50,vmin = 0,vmax = 1000) #vmin = -0.1,vmax = 0.1
    cs.cmap.set_over('#380e0b')
    cs.changed()
    ax6.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    

ani = animation.FuncAnimation(f,func=ani_frame1,frames = za.shape[1],interval = 10,repeat=False)
ani.save('All_combined_Animation_1000mg_wide.mp4',writer=FFwriter)

#%%    
#s[s>100] = 100
baffle_pairs = 6
f = plt.figure(figsize(24,12))
f, (ax3, ax4,ax5,ax6) = plt.subplots(1, 4)
for i in np.arange(0,za.shape[1],1):

    ax3.clear()
    ax3.set_ylim(-0.1,1.5)
    ax3.set_xlim(-0.1,4.1)
    ax3.plot(t[0:i]/60/60/24,I[0:i]/(baffle_pairs*4*A),linewidth = 4)
    ax3.set_ylabel('Current Density (A/$m^2$)',fontsize = 18)
    ax3.set_xlabel('Time (Days)',fontsize = 18)
    ax3.tick_params(axis='both', labelsize=16)
    #ax3.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    ax3.set_title('Average Current Density',fontsize = 20)
    
    ax4.clear()
    ax4.set_ylim(-0.1,3)
    ax4.set_xlim(-0.1,6528)
    ax4.plot(j[...,i],'.')
    ax4.set_ylabel('Current Density (A/$m^2$)',fontsize = 18)
    ax4.set_xlabel('Position',fontsize = 18)
    #ax4.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    ax4.set_title('Positional Current Density',fontsize = 20)
    ax4.tick_params(axis='both', labelsize=16)
    
    ax5.clear()
    ax5.set_ylim(4000,7100)
    ax5.set_xlim(-0.1,6528)
    ax5.plot(za[...,i],'.g')
    ax5.set_ylabel('Biomass (mg/$m^2$)',fontsize = 18)
    ax5.set_xlabel('Position',fontsize = 18)
    #ax5.set_title('Time = {time:.2f} (Days)'.format(time = t[i]/TC_day),fontsize = 20)
    ax5.set_title('Positional Biomass',fontsize = 20)
    ax5.tick_params(axis='both', labelsize=16)
    
    # ax3.clear()
    # ax4.clear()
    # ax5.clear()
    ax6.clear()
    cs = ax6.contourf(xx, yy, s[:, :, i], cmap='coolwarm',levels = 50,vmin = 0,vmax = 1000) #vmin = -0.1,vmax = 0.1
    ax6.set_ylabel('y',fontsize = 18,rotation = 0)
    ax6.set_xlabel('x',fontsize = 18)
    cs.cmap.set_over('#380e0b')
    cs.changed()
    ax6.set_title('Substrate Contour at t = {time:.2f} Days'.format(time = t[i]/TC_day),fontsize = 20)
    plt.tight_layout()
    ax6.tick_params(axis='both', labelsize=16)
    save_figure('hrt_vs_cd_substrate_removal_Y_1000_'+str(i))
    


#%%




for i in np.arange(1,7,1):
    if i == 10:
        data_file_name =  'Y_22_s_500_RT_7_days_HRT_12.npz'
    else:
        data_file_name =  'Y_22_s_500_RT_7_days_HRT_2_baffles_{}.npz'.format(i)
    psi = np.zeros((nx+1,ny+1))  # This is the stream function
    boundary = np.zeros((nx+1,ny+1))  # We set this to 1 on all boundary points
    boundary[0,:] = 1
    boundary[-1,:] = 1
    boundary[:,0] = 1
    boundary[:,-1] = 1
    
    edges = boundary
    psi[0,0:ny+3] = flux
    psi[:,-1] = flux
    psi,boundary,in_out_points,in_start,out_start = system.influent_effluent_regions(int(i), baffle_length, dy*18, psi, boundary, flux)
        
  #  s,za,mox,mred,j = loading_data(data_file_name,['Output','Model_Steady_States','Paper_Results'])
    s,za,mox,mred,j = loading_data(data_file_name,['Output','temp'])
    loaded_baffle_pairs[i - 1] = 2*i
    effluent_substrate[i - 1] = np.mean(s[-2,out_start:out_start+in_out_points],0)
    average_cd[i-1] = j.mean()
    

fig, ax1 = plt.subplots(figsize=(14,10))

color = 'tab:blue'
ax1.set_xlabel('Number of Cartridges',fontsize = 30)
ax1.set_ylabel('Current density (A/m$^2$)', color=color,fontsize = 30)
ax1.plot(loaded_baffle_pairs, average_cd,'.-' ,color=color,linewidth = 5,markersize = 30)
ax1.tick_params(axis='y',labelcolor=color,size = 30)
ax1.tick_params(axis='x',size = 30)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)  

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_xlabel(' Number of Cartridges')
ax2.set_ylabel('Substrate removal %', color=color,fontsize = 30)
ax2.plot(loaded_baffle_pairs, (500-effluent_substrate)/(5),'--.' ,color=color,linewidth = 5,markersize = 30)
ax2.tick_params(axis='y', labelcolor=color,size = 30)
plt.yticks(fontsize = 30)  
plt.title('',fontsize = 30)
plt.tight_layout()
save_figure('hrt_2_cart_comparison_average')

#%%
d_hrt = np.array([0.014842567856412425,0.09988001453202784, 0.24917194967126863, 0.4973472960987408, 1.008407009391439, 2.0465465352345307, 6.111009884629737,18.27142018445394]) 
d_cd = np.array([1.3661396438998314, 1.3145930136929973, 1.2137433693819648, 1.1146336198034459, 1.1008583823348, 1.1698399399646366, 1.0791906533985771, 1.1216497388872897])
d_removal = np.array([0, 8.889550988105498, 24.02364935639301, 51.757826866227504, 58.47275435860432, 54.30674332534159, 75.88017038709529,82.36865994739416])

fig, ax1 = plt.subplots(figsize=(14,18))

color = 'tab:blue'
ax1.set_xlabel('log HRT (Days)',fontsize = 50)
ax1.set_ylabel('Current density (A/m$^2$)', color=color,fontsize = 50)
ax1.plot(d_hrt, d_cd,'.-' ,color='b',linewidth = 5,markersize = 40)
ax1.tick_params(axis='y',labelcolor=color,size = 40)
ax1.tick_params(axis='x',size = 40)
plt.xticks(fontsize = 40)
plt.yticks(fontsize = 40)  

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_xlabel(' log HRT (Days)',fontsize = 50)
ax2.set_ylabel('Substrate removal %', color=color,fontsize = 50)
ax2.semilogx(d_hrt, d_removal,'--.' ,color='r',linewidth = 5,markersize = 40)
ax2.tick_params(axis='y', labelcolor=color,size = 40)
plt.xticks(fontsize = 40)
plt.yticks(fontsize = 40)          
#plt.title('cd and j vs HRT, Mediator yield = 36',fontsize = 30)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()





hrts = [1,2,6,12,24,48,72,144,240]
effluent_substrate = np.zeros(9)
average_cd = np.zeros(9)
loaded_hrt = np.zeros(9)
pos = 0

psi,boundary,in_out_points,in_start,out_start = system.influent_effluent_regions(5, baffle_length, dy*18, psi, boundary, flux)
from pathlib import Path
data_path = Path(os.getcwd(),"Output","Model_Steady_States","Y_40")

for ij in hrts:
    try:
        
        data_file_name =  'Y_22_s_500_RT_7_days_HRT_{}.npz'.format(ij)
        s,za,mox,mred,j = loading_data(data_file_name,['Output','Model_Steady_States','Paper_Results'])
        effluent_substrate[pos] = np.mean(s[-2,out_start:out_start+in_out_points],0)
        average_cd[pos] = j.mean()
        loaded_hrt[pos] = ij
        pos += 1
    
    except:
        import warnings
        print("Data was unable to be loaded for hrt of {}, moving onto next retention time ".format(ij))
        
fig, ax1 = plt.subplots(figsize=(14,18))

color = 'tab:blue'
ax1.set_xlabel('log HRT (Days)',fontsize = 50)
ax1.set_ylabel('Current density (A/m$^2$)', color=color,fontsize = 50)
ax1.semilogx(loaded_hrt/24, average_cd,'.-' ,color=color,linewidth = 5,markersize = 40)
ax1.plot(d_hrt, d_cd,'s-' ,color=color,linewidth = 5,markersize = 20,markerfacecolor='none',markeredgewidth=2)
ax1.tick_params(axis='y',labelcolor=color,size = 40)
ax1.tick_params(axis='x',size = 40)
plt.xticks(fontsize = 40)
plt.yticks(fontsize = 40)  

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_xlabel(' log HRT (Days)',fontsize = 50)
ax2.set_ylabel('Substrate removal %', color=color,fontsize = 50)
ax2.semilogx(loaded_hrt/24, (500-effluent_substrate)/5,'--.' ,color=color,linewidth = 5,markersize = 40)
ax2.semilogx(d_hrt, d_removal,'s--' ,color=color,linewidth = 5,markersize = 20,markerfacecolor='none',markeredgewidth=2)
ax2.tick_params(axis='y', labelcolor=color,size = 40)

  # instantiate a second axes that shares the same x-axis

# color = 'tab:blue'
# ax2.set_ylabel('cd', color=color)  # we already handled the x-label with ax1
# ax2.semilogx(hrt_vec/24, j_vec,'.--' ,color=color,linewidth = 5,markersize = 30)
# ax2.tick_params(axis='y',labelcolor=color)
plt.xticks(fontsize = 40)
plt.yticks(fontsize = 40)          
#plt.title('cd and j vs HRT, Mediator yield = 36',fontsize = 30)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
save_figure('hrt_vs_cd_substrate_removal_experimental')

    


#%%




for _ in hrts:
    try:
        data_storage_path = Path(os.getcwd(),"Output","Model_Steady_States","HRT_{}_Volume_{}_Baffles{}".format(int(_),int(Vol),int(baffle_pairs*2)))
        data_file_name =  'Combined_data_RunTime_7day_'+ anode_side+'.npz'
        loaded_data = np.load(data_storage_path /data_file_name)
        j = loaded_data['cd']
        s = loaded_data['substrate']
        vdata = loaded_data['data']
        average_cd[pos] = np.mean(j)
        effluent_substrate[pos] = np.mean(s[-2,out_start:out_start+in_out_points],0)
        loaded_hrt[pos] = _
        pos +=1
    except:
        import warnings
        print(_)
        warnings.warn("Data was unable to be loaded, moving onto next retention time ")
        
    
data_storage_path = Path(os.getcwd(),"Output","Model_Steady_States","HRT_{}_Volume_{}_Baffles{}".format(int(hrt/TC_hour),int(Vol),int(baffle_pairs*2)))
data_file_name =  'Combined_data_RunTime_7day_'+ anode_side+'.npz'
loaded_data = np.load(data_storage_path /data_file_name)
j = loaded_data['cd']
s = loaded_data['substrate']
vdata = loaded_data['data']

fig, ax1 = plt.subplots(figsize=(14,18))

color = 'tab:blue'
ax1.set_xlabel('log HRT (Days)',fontsize = 30)
ax1.set_ylabel('Current density (A/m$^2$)', color=color,fontsize = 30)
ax1.semilogx(loaded_hrt/24, average_cd,'.--' ,color=color,linewidth = 5,markersize = 30)
ax1.tick_params(axis='y',labelcolor=color,size = 30)
ax1.tick_params(axis='x',size = 30)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)  

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_xlabel(' log HRT (Days)',fontsize = 30)
ax2.set_ylabel('Substrate removal %', color=color,fontsize = 30)
ax2.semilogx(loaded_hrt/24, (500-effluent_substrate)/5,'--.' ,color=color,linewidth = 5,markersize = 30)
ax2.tick_params(axis='y', labelcolor=color,size = 30)

  # instantiate a second axes that shares the same x-axis

# color = 'tab:blue'
# ax2.set_ylabel('cd', color=color)  # we already handled the x-label with ax1
# ax2.semilogx(hrt_vec/24, j_vec,'.--' ,color=color,linewidth = 5,markersize = 30)
# ax2.tick_params(axis='y',labelcolor=color)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)          
#plt.title(title,fontsize = 26)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
#save_figure('hrt_vs_cd_substrate_removal')


# plt.figure(figsize = (10,10))
# plt.loglog(hrt_vec/24,i_totals/7,'o--g',linewidth = 5,markersize =20)
# plt.xlabel('Log HRT (Days)',fontsize = 20)
# plt.ylabel('Average Daily Current A/D',fontsize = 20)
# plt.xticks(fontsize = 20)
# plt.yticks(fontsize = 20)
# plt.title('Average Daily current vs HRT',fontsize = 26)


# plt.figure()
# plt.subplot(211)
# plot_positional_data(x,j,bio_loc,side = 'Right')
# plt.subplot(212)
# plot_contour(xx,yy,s,ux_data = 0,uy_data = 0,scatter_positions = bio_loc,scatter_data = vdata[0],title=anode_side )

# plt.figure(figsize = (10,12))
# plot_contour(xx,yy,s,ux_data = 0,uy_data = 0,scatter_positions = bio_loc,scatter_data = vdata[0],title='HRT = {}'.format(int(hrt/TC_hour)) )


# data_storage_path = Path(os.getcwd(),"Output","Model_Steady_States","HRT_{}_Volume_{}_Baffles{}".format(int(48),int(Vol),int(baffle_pairs*2)))
# data_file_name =  'Combined_data_RunTime_7day_'+ anode_side+'.npz'
# loaded_data = np.load(data_storage_path /data_file_name)
# s = loaded_data['substrate']
# vdata = loaded_data['data']
# plt.figure()
# plot_contour(xx,yy,s,new_fig = True,scatter_positions = bio_loc,scatter_data = vdata[0],title = 'HRT = 48 hours')

#%%
# new_vdata = np.zeros((6,1088))
# for i in np.arange(1,7):
#     print(i)
#     nx = 300
#     ny = 300
#     nz = 1
#     # Size
#     Lx = 0.32 # length in m
#     Ly = 0.45 # length in m
#     Lz = 0.25 # length in m
#     # HRT
#     hrt = 6  # hrt set in any given time units TC is the conversion to seconds
#     hrt *= TC_hour  # The hydraulic retention time converted into seconds
#     #Baffle param
#     baffle_length = 91/100 ##### This is the fraction of the tank a baffle takes up in x
#     baffle_pairs = int(i) #### Number of baffle pairs (RHS+LHS) = 1 pair.
#     # Runtime for reactor system
#     RT  = 7
#     RT *= TC_day
#     dt_max = 4
#     k = 20 # Store data every k min, 1 is every min, 2 every 2 min
#     D = 1/hrt
#     stream = False
#     anode_side = 'all' # can either be 'all', 'influent', 'effluent'
#     Vol = int(Lx*Ly*Lz*1e3)
    
    
#     dx = Lx/nx
#     dy = Ly/ny
#     x = np.linspace(0,Lx,nx).T
#     y = np.linspace(0,Ly,ny).T
#     [yy,xx] = np.meshgrid(np.linspace(dy/2,Ly-dy/2,ny),np.linspace(dx/2,Lx-dx/2,nx))
    
#     system = domain(Lx, Ly, Lz,nx = nx,ny = ny)
    
    
#     flux = (Lx*Ly)/hrt # The "area flux" through the system
#     nxy = nx*ny
#     nxy_one = (nx+1)*(ny+1)
#     psi = np.zeros((nx+1,ny+1))  # This is the stream function
#     boundary = np.zeros((nx+1,ny+1))  # We set this to 1 on all boundary points
#     boundary[0,:] = 1
#     boundary[-1,:] = 1
#     boundary[:,0] = 1
#     boundary[:,-1] = 1
    
#     edges = boundary
#     psi[0,0:ny+3] = flux
#     psi[:,-1] = flux
#     psi,boundary,in_out_points,in_start,out_start = system.influent_effluent_regions(baffle_pairs, baffle_length, dy*18, psi, boundary, flux)
#     bdata = psi[boundary == 1]
#     external = np.zeros(boundary.shape)
#     external[0,:] = 1
#     external[-1,:] = 1
#     external[:,0] = 1
#     external[:,-1] = 1
#     internal = boundary - external
#     bio_loc = np.zeros(boundary.shape)
#     bio_loc[:,0:ny] += internal[:,1:ny+1]
#     bio_loc[:,1:ny+1] += internal[:,0:ny]
#     bio_loc = bio_loc[0:nx,0:ny]
#     if baffle_length==0 and baffle_pairs == 0:
#         bio_loc[1:-1,-2] = 1
#         bio_loc[1:-1,1] = 1
#     file_a = 'hrt' + str(hrt).replace('.','_') + '_nx' +str(nx) + '_ny' +str(ny)+ '_Lx' +str(Lx) + '_Ly' +str(Ly)+'_pairs'+str(baffle_pairs)+'_width'+str(np.round(baffle_length,decimals=1)).replace('.','_')+'.csv'
#     file_x = 'Ux_' + file_a
#     file_y = 'Uy_' + file_a
#     data_folder = Path(os.getcwd(),"Output","Velocity")
#     ux = np.genfromtxt(data_folder / file_x, delimiter=',')
#     uy = np.genfromtxt(data_folder / file_y, delimiter=',')
    
    
#     positional = np.nonzero(np.mean(bio_loc,0))
#     switch = np.zeros(bio_loc.shape)
#     if anode_side == 'influent':
#         switch[:,positional[0][0:20:2]] = 1 #:2 to alternate front and back on off 1:20:2 is back, 0:20:2 is front
#         bio_loc *= switch
#         print('influent facing anodes are active')
#     elif anode_side == 'effluent':
#         switch[:,positional[0][1:20:2]] = 1 #:2 to alternate front and back on off 1:20:2 is back, 0:20:2 is front
#         bio_loc *= switch
#         print('effluent facing anodes are active')
#     else:
#         anode_side = 'all'
#         print('All anodes are active')
#     single_switch = np.zeros(bio_loc.shape)
#     single_positional = np.nonzero(np.mean(bio_loc,0))
#     single_switch[:,single_positional[0][0:4]] = 1
#     single_data = bio_loc *single_switch
#     pos = baffle_pairs - 1
#     data_storage_path = Path(os.getcwd(),"Output","Model_Steady_States","HRT_{}_Volume_{}_Baffles{}".format(int(6),int(Vol),int(baffle_pairs*2)))
#     data_file_name =  'Combined_data_RunTime_7day_'+ anode_side+'.npz'
#     loaded_data = np.load(data_storage_path /data_file_name)
#     j = loaded_data['cd']
#     s = loaded_data['substrate']
#     vdata = loaded_data['data']
#     full_data = np.zeros(bio_loc.shape)
#     full_data[bio_loc==1] = vdata[0]
#     single_data *= full_data
#     new_vdata[int(i-1),:] = single_data[single_data != 0]
#     average_cd[pos] = np.mean(j)
#     effluent_substrate[pos] = np.mean(s[-2,out_start:out_start+in_out_points],0)
#     loaded_hrt[pos] = baffle_pairs


nx = 300
ny = 300
nz = 1
# Size
Lx = 0.32 # length in m
Ly = 0.45 # length in m
Lz = 0.25 # length in m
# HRT
hrt = 12 # hrt set in any given time units TC is the conversion to seconds
hrt *= TC_hour  # The hydraulic retention time converted into seconds
#Baffle param
baffle_length = 91/100 ##### This is the fraction of the tank a baffle takes up in x #### Number of baffle pairs (RHS+LHS) = 1 pair.
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
plt.figure(figsize = (14,10))
for ij in np.arange(1,5,1): #7
    baffle_pairs = int(ij)
    early_linesty = ['b-','r-','g-','k-']
    late_linesty = ['b--','r--','g--','k--']
    
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
    psi,boundary,in_out_points,in_start,out_start = system.influent_effluent_regions(baffle_pairs, baffle_length, dy*18, psi, boundary, flux)
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
    ux = np.genfromtxt(data_folder / file_x, delimiter=',')
    uy = np.genfromtxt(data_folder / file_y, delimiter=',')
    
    
    positional = np.nonzero(np.mean(bio_loc,0))
    switch = np.zeros(bio_loc.shape)    
    pos = baffle_pairs - 1
    data_storage_path = Path(os.getcwd(),"Output","Model_Steady_States","HRT_{}_Volume_{}_Baffles{}".format(int(6),int(Vol),int(baffle_pairs*2)))
    data_file_name =  'Combined_data_RunTime_7day_'+ anode_side+'.npz'
    loaded_data = np.load(data_storage_path /data_file_name)
    j = loaded_data['cd']
    s = loaded_data['substrate']
    vdata = loaded_data['data']
    if baffle_pairs == 50:
        data_file_name =  'Y_22_s_500_RT_7_days_HRT_12.npz'
    else:
        data_file_name =  'Y_22_s_500_RT_7_days_HRT_72_baffles_{}.npz'.format(baffle_pairs)
    #s,z,mox,mred,j = loading_data(data_file_name,['Output','Model_Steady_States','Paper_Results'])
    s,z,mox,mred,j = loading_data(data_file_name,['Output','temp'])
    
    reduced_bio_loc = np.zeros(bio_loc.shape)
    first_array = np.arange(0,4,1)
    second_array = np.arange(-4,0,1)
    positions_list = list(first_array) + list(second_array)
    switch[:,positional[0][positions_list]] = 1 #[0:4:1]
    reduced_bio_loc = bio_loc*switch
    distr = np.zeros(bio_loc.shape)
    distr[bio_loc == 1] = j
    distr *= reduced_bio_loc
    j = distr[reduced_bio_loc == 1]
    
    switch = np.zeros(bio_loc.shape)    
    switch[:,positional[0][0]] = 1
    bio_loc_front_left = bio_loc*switch 
    j_fl = distr[bio_loc_front_left == 1]
    switch = np.zeros(bio_loc.shape)    
    switch[:,positional[0][1]] = 1
    bio_loc_back_left = bio_loc*switch 
    j_bl = distr[bio_loc_back_left == 1]
    switch = np.zeros(bio_loc.shape)    
    switch[:,positional[0][2]] = 1
    bio_loc_front_right = bio_loc*switch 
    switch = np.zeros(bio_loc.shape)   
    j_fr = distr[bio_loc_front_right == 1]
    switch[:,positional[0][3]] = 1
    bio_loc_back_right = bio_loc*switch 
    j_br = distr[bio_loc_back_right == 1]
    
    switch = np.zeros(bio_loc.shape)    
    switch[:,positional[0][-4]] = 1
    bio_loc_front_left = bio_loc*switch 
    j2_fl = distr[bio_loc_front_left == 1]
    switch = np.zeros(bio_loc.shape)    
    switch[:,positional[0][-3]] = 1
    bio_loc_back_left = bio_loc*switch 
    j2_bl = distr[bio_loc_back_left == 1]
    switch = np.zeros(bio_loc.shape)    
    switch[:,positional[0][-2]] = 1
    bio_loc_front_right = bio_loc*switch 
    switch = np.zeros(bio_loc.shape)   
    j2_fr = distr[bio_loc_front_right == 1]
    switch[:,positional[0][-1]] = 1
    bio_loc_back_right = bio_loc*switch 
    j2_br = distr[bio_loc_back_right == 1]
    
    

    plt.subplot(221)
    plot_positional_data(x,j_bl,bio_loc_back_left,side = 'Left',legend_on= False,xlab = '',ylab = 'Current Density (A/m$^2$)',title = 'Effluent Left',sty = early_linesty[ij-1])
    plot_positional_data(x,j2_bl,bio_loc_back_left,side = 'Left',legend_on= False,xlab = '',sty = late_linesty[ij-1],ylab = 'Current Density (A/m$^2$)',title = 'Effluent Left')

    plt.xticks(np.arange(min(x), max(x), 0.1))
    if baffle_pairs == 4:
        plt.ylim([0,1.5])
        plt.legend(['2','2','4','4','6','6','8','8'],fontsize = 16, ncol=4)
    plt.subplot(222)
    plot_positional_data(x,j_br,bio_loc_back_right,side = 'Right',legend_on= False,xlab = '',title = 'Effluent Right',sty = early_linesty[ij-1])
    plot_positional_data(x,j2_br,bio_loc_back_right,side = 'Right',legend_on= False,xlab = '',sty = late_linesty[ij-1],title = 'Effluent Right')
    plt.xticks(np.arange(min(x), max(x), 0.1))
    if baffle_pairs == 4:
        plt.ylim([0,1.5])
        plt.legend(['2','2','4','4','6','6','8','8'],fontsize = 16, ncol=4)
        #plt.legend(['2','4','6','8','10','12'],fontsize = 16, ncol=3)
    plt.subplot(223)
    plot_positional_data(x,j_fl,bio_loc_front_left,side = 'Left',legend_on= False,xlab = 'x',ylab = 'Current Density (A/m$^2$)',title = 'Influent Left',sty = early_linesty[ij-1])
    plot_positional_data(x,j2_fl,bio_loc_front_left,side = 'Left',legend_on= False,xlab = 'x',sty = late_linesty[ij-1],ylab = 'Current Density (A/m$^2$)',title = 'Influent Left')
    plt.xticks(np.arange(min(x), max(x), 0.1))
    if baffle_pairs == 4:
        plt.ylim([0,1.5])
        plt.legend(['2','2','4','4','6','6','8','8'],fontsize = 16, ncol=4)
        #plt.legend(['2','4','6','8','10','12'],fontsize = 16, ncol=3)
    plt.subplot(224)
    plot_positional_data(x,j_fr,bio_loc_front_right,side = 'Right',legend_on= False,xlab = 'x',title = 'Influent Right',sty = early_linesty[ij-1])
    plot_positional_data(x,j2_fr,bio_loc_front_right,side = 'Right',legend_on= False,xlab = 'x',sty = late_linesty[ij-1],title = 'Influent Right')
    plt.xticks(np.arange(min(x), max(x), 0.1))
    if baffle_pairs == 4:
        plt.ylim([0,1.5])
        plt.legend(['2','2','4','4','6','6','8','8'],fontsize = 16, ncol=4)
        #plt.legend(['2','4','6','8','10','12'],fontsize = 16, ncol=3,loc =3)
    back_mean = (j2_bl.mean() + j2_br.mean() + j_bl.mean() + j_br.mean())/4
    front_mean = (j2_fl.mean() + j2_fr.mean() + j_fl.mean() + j_fr.mean())/4
    
    print('mean values for {} baffles are \n {} Effluent facing \n {} Influent Facing'.format(2*baffle_pairs,back_mean, front_mean))
    
plt.tight_layout()
save_figure('Anode_sides_current_density_reshaped_first_last_anodes')

#%%
def anode_surface_sum(value, location):
    location = np.array(location, dtype=bool)
    temp = np.zeros(location.shape)
    temp[location] = value
    anode_surface_sum = np.sum(temp,0)
    return anode_surface_sum


def anode_surface_sum_repeated(value, location, additional_pref = 1):
    location = np.array(location, dtype=bool)
    temp = anode_surface_sum(value, location)

    # location = np.array(location, dtype=bool)
    # temp = np.zeros(location.shape)
    # temp[location] = value

    anode_total = np.array([additional_pref*temp, ] * nx)
    anode_total = anode_total[location]
    return anode_total

s,za,mox,mred,j = loading_data('Y_40_s_500_RT_7_days_HRT_24.npz',['Output','Model_Steady_States','Paper_Results'])
Ai = 0.0002666666666666667
A = 0.0728
I_anode_repeated = anode_surface_sum_repeated(j*0.0002666666666666667,bio_loc)

eta_act = 1.4*(26*0.0728) - R*T/(2*F)*np.log((mox+mred)/mred) - I_anode_repeated * 26
print(eta_act.mean())

#%%
hrt_d = np.array([0.09988001453202784, 0.24917194967126863, 0.4973472960987408, 1.008407009391439, 2.0465465352345307, 6.111009884629737,18.27142018445394]) 
cd_d = np.array([ 1.3145930136929973, 1.2137433693819648, 1.1146336198034459, 1.1008583823348, 1.1698399399646366, 1.0791906533985771, 1.1216497388872897])
removal_d = np.array([8.889550988105498, 24.02364935639301, 51.757826866227504, 58.47275435860432, 54.30674332534159, 75.88017038709529,82.36865994739416])
plt.plot(hrt,removal)
plt.plot(hrt,cd)

fig, ax1 = plt.subplots(figsize=(14,18))


hrts = [1,2,6,12,24,48,72,144,240]
effluent_substrate = np.zeros(9)
average_cd = np.zeros(9)
loaded_hrt = np.zeros(9)
pos = 0

psi,boundary,in_out_points,in_start,out_start = system.influent_effluent_regions(5, baffle_length, dy*18, psi, boundary, flux)
from pathlib import Path
data_path = Path(os.getcwd(),"Output","Model_Steady_States","Y_40")

for ij in hrts:
    try:
        
        data_file_name =  'Y_22_s_500_RT_7_days_HRT_{}.npz'.format(ij)
        s,za,mox,mred,j = loading_data(data_file_name,['Output','Model_Steady_States','Paper_Results'])
        effluent_substrate[pos] = np.mean(s[-2,out_start:out_start+in_out_points],0)
        average_cd[pos] = j.mean()
        loaded_hrt[pos] = ij
        pos += 1
    
    except:
        import warnings
        print("Data was unable to be loaded for hrt of {}, moving onto next retention time ".format(ij))

color = 'tab:blue'
ax1.set_xlabel('log HRT (Days)',fontsize = 50)
ax1.set_ylabel('Current density (A/m$^2$)', color=color,fontsize = 50)
ax1.semilogx(loaded_hrt/24, average_cd,'.-' ,color=color,linewidth = 5,markersize = 40)
ax1.plot(hrt_d, cd_d,'.-' ,color='g',linewidth = 5,markersize = 40)
ax1.tick_params(axis='y',labelcolor=color,size = 40)
ax1.tick_params(axis='x',size = 40)
plt.xticks(fontsize = 40)
plt.yticks(fontsize = 40)  

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_xlabel(' log HRT (Days)',fontsize = 50)
ax2.set_ylabel('Substrate removal %', color=color,fontsize = 50)
ax2.semilogx(loaded_hrt/24, (500-effluent_substrate)/5,'--.' ,color=color,linewidth = 5,markersize = 40)
ax2.semilogx(hrt_d, removal_d,'--.' ,color='g',linewidth = 5,markersize = 40)
ax2.tick_params(axis='y', labelcolor=color,size = 40)

  # instantiate a second axes that shares the same x-axis

# color = 'tab:blue'
# ax2.set_ylabel('cd', color=color)  # we already handled the x-label with ax1
# ax2.semilogx(hrt_vec/24, j_vec,'.--' ,color=color,linewidth = 5,markersize = 30)
# ax2.tick_params(axis='y',labelcolor=color)
plt.xticks(fontsize = 40)
plt.yticks(fontsize = 40)          
#plt.title('cd and j vs HRT, Mediator yield = 36',fontsize = 30)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
#save_figure('hrt_vs_cd_substrate_removal_Y_22')

# color = 'tab:blue'
# ax1.set_xlabel('log HRT (Days)',fontsize = 50)
# ax1.set_ylabel('Current density (A/m$^2$)', color=color,fontsize = 50)
# ax1.plot(hrt_d, cd_d,'.-' ,color='g',linewidth = 5,markersize = 40)
# ax1.tick_params(axis='y',labelcolor=color,size = 40)
# ax1.tick_params(axis='x',size = 40)
# plt.xticks(fontsize = 40)
# plt.yticks(fontsize = 40)  

# ax2 = ax1.twinx()
# color = 'tab:red'
# ax2.set_xlabel(' log HRT (Days)',fontsize = 50)
# ax2.set_ylabel('Substrate removal %', color=color,fontsize = 50)
# ax2.semilogx(hrt_d, removal_d,'--.' ,color='g',linewidth = 5,markersize = 40)
# ax2.tick_params(axis='y', labelcolor=color,size = 40)
# plt.xticks(fontsize = 40)
# plt.yticks(fontsize = 40)          
# #plt.title('cd and j vs HRT, Mediator yield = 36',fontsize = 30)
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# plt.show()




