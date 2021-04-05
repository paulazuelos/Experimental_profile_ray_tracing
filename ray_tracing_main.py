from libraries_import import *
from lens_profile_interpolation import *
from emitter_positionning import *
from ray_tracing_engine import *
from plotting_functions import *
from general_functions import *
# P.Azuelos Idemia R&T team 2020
# Forward ray tracing algorithm;  LQA project 
# study of the effect of the lens shape
# study of the effect of the engraving and color matrix position
# render above lenses with an angle
# compare sperical and parabolic lens profile
#%% input 
    
####### input parameters #####################################################

wavelength=0.6 # wavelength in µm   
n_PC=n_PC_func(wavelength) # refractive index of the polycarbonate

n_air=1 # refractive index of the air
n=n_air/n_PC

### input variable for the theoretical lens profile (choose_data=0)

lens_period=90 # pitch between two lenses
lens_heigth=17 # heigth of the lens (cylindrical profile)
nb_lens=10 # number of lens to be computed

### input variable for the experimental lens profile (choose_data=1)

expected_width=140 # used for smoothing and peak finding of the raw experimental datas
split_lens=math.floor(expected_width/4) # discretization steps per lens check the results for speed optimisation/ not to high to avoid simulating noise
filter_detect_peak=0.2 # power of the filter for the peak detection

### Directory and file name managment

file_path='C40P1 HRX_MM  - Extracted profile.txt'
x_crop_min=100 #x start for the profile calculation
x_crop_max=15000 #x stop for the profile calculation

### variables to drive the software parts

choose_data=0# 0 generate a theoretical lens array, 1 interpolate an experimental lens array, other number don't calculate again
lamination_plate=0 # sometimes we want to check the difference between plates profile and actual laminated lens (it just flip the data along y-axis)
calculation=1 # calculate only when the lens profile and the calculation points are set up
plotting_experimental=0
ray_tracing=0

### Angle of calculation

theta_lim=2*math.asin(n_air/n_PC)*180/math.pi
out_angle_lim=90 # range of angles for the plot
split_out_angle=2*out_angle_lim # for the discretization of angle range
range_angle=np.linspace(-out_angle_lim,out_angle_lim,split_out_angle)

### PSF calculation

PSF=1
dist_center=0 # position the PSF emiter at a certain distance from the center of the lens (e.g: 10 => 10 µm from the lens center axis)
nb_step=25 # discretization of the PSF
y_top_PSF=100
y_bottom_PSF=300

z_depth=np.linspace(y_top_PSF,y_bottom_PSF,nb_step)
Z_scan_PSF=np.zeros((range_angle.size,nb_step))



#%% Calculation and interpolation of the lens profiles

if choose_data==0:
    x_grid, y_grid,F,Y_lens,peakind,nb_lens_real,lens_period_mean=Interp_profile_theo(lens_period,lens_heigth,n_PC,nb_lens,split_lens)
if choose_data==1:
    x_grid, y_grid, F, Y_lens, peakind,nb_lens_real,lens_period_mean=Interp_profile_exp(file_path,x_crop_min,x_crop_max,expected_width,filter_detect_peak,lamination_plate,split_lens,n_PC)

plot_profile(x_grid,y_grid)
#%% Experimental plot of the lens profiles
#if plotting_experimental==1:
X_lens=x_grid    
nb_points_max=100000
nb_cut=8
### input variables for calculation
nb_ray_calc_max=800000
### input variables for plot
colors = np.array([(1, 1, 0), (0, 1, 1), (1, 0, 1),(0.5,0.5,0.5)])




if ray_tracing==1:
    step_calc=5
if ray_tracing==0:
    step_calc=1

for idx, nb_step_calculation_points in enumerate(z_depth):
    ### input variables for the positionning and discretization of the calculation points

    print("step: {0}".format(idx))
    y_lim_calculation=-nb_step_calculation_points
    y_lim_top=math.floor(nb_step_calculation_points/4)
    y_grid=np.linspace(0,y_lim_calculation,nb_step_calculation_points+1)
    inc_x=abs(x_grid[0]-x_grid[1])
    inc_y=abs(y_grid[0]-y_grid[1])


    #############################################################################
    #############################################################################
    #%% Definition of the diffusion points (matrix or engraving position)
    # the points can be positionned anywhere below the lens 
    
    color_def=Emitters_positionning(y_lim_calculation,nb_step_calculation_points,x_grid,y_grid,lens_period,nb_cut,peakind,nb_points_max,PSF,dist_center)
    
    color_def=color_def[~np.all(color_def == 0, axis=1)]
    color_matrix=np.array(color_def[:,0].astype(int))
    colors = np.array([(1, 1, 0), (0, 1, 1), (1, 0, 1),(0.5,0.5,0.5)])
    f2 = plt.figure()
    Font_SIZE = 15
    plt.rc('font', size=Font_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=Font_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=Font_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=Font_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=Font_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=Font_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=Font_SIZE)  # fontsize of the figure title
    
    plt.xlabel('x (um)')
    plt.ylabel('y (um)')
    plt.plot(X_lens, Y_lens, '.')
    x_vect=np.ones(len(y_grid))
   
    plt.scatter(color_def[:,3], color_def[:,4], s=1, c=colors[color_matrix])
    
    x1,x2,y1,y2 = plt.axis()
    #plt.axis((T_grid_min,T_grid_max,y1,90))
    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
    #f2.savefig('measured_conf_{0}.png'.format(name_file), dpi=800, bbox_inches='tight', pad_inches=0)
    plt.show()
    #%% Calculation (forward ray tracing vectors calculation)
    
    if calculation==1:
        ##### Calculate all the forward rays angles direction
        out_rays=Ray_tracing_engine(F,theta_lim,inc_x,nb_ray_calc_max,color_def,X_lens,Y_lens,step_calc,n_PC,n_air,y_lim_top)
        # Estimated the mean color depending on the angle of view and x-position
                    
    #%% Sort all the vectors depending on their angle and x-position 
        
        range_x=np.linspace(np.min(X_lens),np.max(X_lens),nb_lens)
        
        delta_x=range_x[1]-range_x[0]
        delta_angle=range_angle[1]-range_angle[0]
        
        
        value_mean=np.zeros((range_x.size,range_angle.size,4))
        for i in range(0,range_x.size):
            for j in range(0,range_angle.size):
                test=np.where(np.logical_and(out_rays[:,2]>range_angle[j-1], out_rays[:,2]<=range_angle[j]),0,1)
                test+=np.where(np.logical_and(out_rays[:,3]>range_x[i-1], out_rays[:,3]<=range_x[i]),0,1)
                select=out_rays[test==0]
                if select.size:
        
                    
                    T1=select[np.where(select[:,4]==0,0,1)==0,7]
                    T2=select[np.where(select[:,4]==1,0,1)==0,7]
                    T3=select[np.where(select[:,4]==2,0,1)==0,7]
                    T4=select[np.where(select[:,4]==3,0,1)==0,7]
                    
        
                    value_mean[i-1,j-1,0]=np.sum(T1)/delta_angle # first dimension (x (µm)) second dimension (angle (°)) color 0
                    value_mean[i-1,j-1,1]=np.sum(T2)/delta_angle # first dimension (x (µm)) second dimension (angle (°)) color 1
                    value_mean[i-1,j-1,2]=np.sum(T3)/delta_angle # first dimension (x (µm)) second dimension (angle (°)) color 2
                    value_mean[i-1,j-1,3]=np.sum(T4)/delta_angle # first dimension (x (µm)) second dimension (angle (°)) color 3
        
    #%% Plot the results 
    # I must check the intensity value (for sure divide by the number of lenses used to generate the plot)
    # the discretization variables are not the sames for theoretical and experimental lens profiles (can be homogeneized in the futur)
    
        if choose_data==1 and plotting_experimental==1: 

            f1 = plt.figure()
            Font_SIZE = 15
            plt.rc('font', size=Font_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=Font_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=Font_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=Font_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=Font_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=Font_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=Font_SIZE)  # fontsize of the figure title
            
            plt.xlabel('Angle (°)')
            plt.ylabel('Normalized Point Spread Function')
            plt.plot(range_angle, np.sum(value_mean[1:nb_lens_real-2,:,0],axis=0)/max(np.sum(np.sum(value_mean[1:nb_lens_real-2,:,0],axis=0)),1), '-',c=colors[0])
            plt.plot(range_angle, np.sum(value_mean[1:nb_lens_real-2,:,1],axis=0)/max(np.sum(np.sum(value_mean[1:nb_lens_real-2,:,1],axis=0)),1), '-',c=colors[1])
            plt.plot(range_angle, np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)/max(np.sum(np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)),1), '-',c=colors[2])
            plt.plot(range_angle, np.sum(value_mean[1:nb_lens_real-2,:,3],axis=0)/max(np.sum(np.sum(value_mean[1:nb_lens_real-2,:,3],axis=0)),1), '-',c=colors[3])

            f1.savefig('x_theta_render_measured_{0}_split_lens_{1}.png'.format(name_file,split_lens).format(name_file), dpi=800, bbox_inches='tight', pad_inches=0) # For illustration
            color_matrix=np.array(color_def[:,0].astype(int))
            
            X=X_lens[hollow_ind[1]:hollow_ind[2], np.newaxis]
            Y=Y_lens[hollow_ind[1]:hollow_ind[2], np.newaxis]
            
            A = np.hstack([X**2, X * Y, Y**2, X, Y])
            b = np.ones_like(X)
            x = np.linalg.lstsq(A, b)[0].squeeze()
            
            x_coord = np.linspace(X_lens[hollow_ind[1]],X_lens[hollow_ind[2]],100)
            y_coord = np.linspace(np.amin(Y_lens[hollow_ind[1]:hollow_ind[2]]),np.amax(Y_lens[hollow_ind[1]:hollow_ind[2]]),100)
            X_coord, Y_coord = np.meshgrid(x_coord, y_coord)
            Z_coord = x[0] * X_coord ** 2 + x[1] * X_coord * Y_coord + x[2] * Y_coord**2 + x[3] * X_coord + x[4] * Y_coord
            
            if ray_tracing==1:
                
                color_def=color_def[~np.all(color_def == 0, axis=1)]
                color_matrix=np.array(color_def[:,0].astype(int))
                colors = np.array([(1, 1, 0), (0, 1, 1), (1, 0, 1),(0.5,0.5,0.5)])
                f4 = plt.figure()
                Font_SIZE = 15
                
                plt.rc('font', size=Font_SIZE)          # controls default text sizes
                plt.rc('axes', titlesize=Font_SIZE)     # fontsize of the axes title
                plt.rc('axes', labelsize=Font_SIZE)    # fontsize of the x and y labels
                plt.rc('xtick', labelsize=Font_SIZE)    # fontsize of the tick labels
                plt.rc('ytick', labelsize=Font_SIZE)    # fontsize of the tick labels
                plt.rc('legend', fontsize=Font_SIZE)    # legend fontsize
                plt.rc('figure', titlesize=Font_SIZE)  # fontsize of the figure title
                
                plt.xlabel('x (um)')
                plt.ylabel('y (um)')
                plt.plot(X_lens, Y_lens, '.')
                x_vect=np.ones(len(y_grid))
    
                plt.scatter(color_def[:,3], color_def[:,4], s=1, c=colors[color_matrix])
                theta_one_lens=0.8*(math.pi/2-math.atan(abs(y_lim_calculation)/(lens_period_mean/2)))*(180/math.pi)
                for i in range(0,out_rays.shape[0]):
                    if out_rays[i,8]<X_lens[-1] and out_rays[i,8]>X_lens[0] and abs(out_rays[i,0])<theta_one_lens:
                        
                        plt.plot([out_rays[i,3],out_rays[i,5]], [out_rays[i,9],out_rays[i,10]],'-b', linewidth=0.25)
                        plt.plot([out_rays[i,5],out_rays[i,8]], [out_rays[i,10],out_rays[i,11]],'-b', linewidth=0.25)
                    
                x1,x2,y1,y2 = plt.axis()
                #plt.axis((T_grid_min,T_grid_max,y1,90))
                plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
                f4.savefig('measured_conf_ray_tracing_{0}_split_lens_{1}.png'.format(name_file,split_lens), dpi=800, bbox_inches='tight', pad_inches=0)
                plt.show()
            
            
        if choose_data==0 and plotting_experimental==1:
            
            f1 = plt.figure()
            Font_SIZE = 15
            plt.rc('font', size=Font_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=Font_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=Font_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=Font_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=Font_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=Font_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=Font_SIZE)  # fontsize of the figure title
            
            plt.xlabel('Angle (°)')
            plt.ylabel('Normalized Point Spread Function')
            plt.plot(range_angle, np.sum(value_mean[1:nb_lens-2,:,0],axis=0)/np.sum(np.sum(value_mean[1:nb_lens-2,:,0],axis=0)), '-',c=colors[0])
            plt.plot(range_angle, np.sum(value_mean[1:nb_lens-2,:,1],axis=0)/np.sum(np.sum(value_mean[1:nb_lens-2,:,1],axis=0)), '-',c=colors[1])
            plt.plot(range_angle, np.sum(value_mean[1:nb_lens-2,:,2],axis=0)/np.sum(np.sum(value_mean[1:nb_lens-2,:,2],axis=0)), '-',c=colors[2])
            plt.plot(range_angle, np.sum(value_mean[1:nb_lens-2,:,3],axis=0)/np.sum(np.sum(value_mean[1:nb_lens-2,:,3],axis=0)), '-',c=colors[3])
            f1.savefig('x_theta_render_measured_spherical_profil_split_lens_{0}.png'.format(split_lens), dpi=800, bbox_inches='tight', pad_inches=0) # For illustration
            color_matrix=np.array(color_def[:,0].astype(int))
            
            if ray_tracing==1:
            
                color_def=color_def[~np.all(color_def == 0, axis=1)]
                color_matrix=np.array(color_def[:,0].astype(int))
                colors = np.array([(1, 1, 0), (0, 1, 1), (1, 0, 1),(0.5,0.5,0.5)])
                f4 = plt.figure()
                Font_SIZE = 15
                
                plt.rc('font', size=Font_SIZE)          # controls default text sizes
                plt.rc('axes', titlesize=Font_SIZE)     # fontsize of the axes title
                plt.rc('axes', labelsize=Font_SIZE)    # fontsize of the x and y labels
                plt.rc('xtick', labelsize=Font_SIZE)    # fontsize of the tick labels
                plt.rc('ytick', labelsize=Font_SIZE)    # fontsize of the tick labels
                plt.rc('legend', fontsize=Font_SIZE)    # legend fontsize
                plt.rc('figure', titlesize=Font_SIZE)  # fontsize of the figure title
                
                plt.xlabel('x (um)')
                plt.ylabel('y (um)')
                plt.plot(X_lens, Y_lens, '.')
                x_vect=np.ones(len(y_grid))
        
                plt.scatter(color_def[:,3], color_def[:,4], s=1, c=colors[color_matrix])
                theta_one_lens=0.8*(math.pi/2-math.atan(abs(y_lim_calculation)/(lens_period_mean/2)))*(180/math.pi)
                for i in range(0,out_rays.shape[0]):
                    if out_rays[i,8]<X_lens[-1] and out_rays[i,8]>X_lens[0] and abs(out_rays[i,0])<theta_one_lens:
                        
                        plt.plot([out_rays[i,3],out_rays[i,5]], [out_rays[i,9],out_rays[i,10]],'-b', linewidth=0.25)
                        plt.plot([out_rays[i,5],out_rays[i,8]], [out_rays[i,10],out_rays[i,11]],'-b', linewidth=0.25)
                    
                x1,x2,y1,y2 = plt.axis()
                #plt.axis((T_grid_min,T_grid_max,y1,90))
                plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
                f4.savefig('measured_conf_ray_tracing_spherical_profil_split_lens_{0}.png'.format(split_lens), dpi=800, bbox_inches='tight', pad_inches=0)
                plt.show()
                
        if PSF==0:
            
            mean_2=np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)
            mean_1=np.sum(value_mean[1:nb_lens_real-2,:,1],axis=0)
            contrast=(mean_2-mean_1)/(mean_2+mean_1)
            
            Z_scan_PSF[:,idx]=contrast
        
        if PSF==1:
            Z_scan_PSF[:,idx]=np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)/max(np.sum(np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)),1)
#%%            
z_depth_map,range_angle_map=np.meshgrid(z_depth,range_angle)

#%% PSF Map plot

fig5, ax2 = plt.subplots(constrained_layout=True)
origin = 'lower'

CS = ax2.contourf(range_angle_map, z_depth_map,Z_scan_PSF, 100, cmap=plt.cm.bwr, origin=origin)

if PSF==1:
    focal_pos_z=z_depth_map[np.unravel_index(np.argmax(Z_scan_PSF, axis=None), Z_scan_PSF.shape)]
    focal_pos_theta=range_angle_map[np.unravel_index(np.argmax(Z_scan_PSF, axis=None), Z_scan_PSF.shape)]
    ax2.scatter(focal_pos_theta, focal_pos_z,s=500,marker="x", color="black")
    print("focal pos : theta={0} z={1}".format(focal_pos_theta,focal_pos_z))
    ax2.set_title('PSF with z-position')
    ax2.set_xlabel('Angle(°)')
    ax2.set_ylabel('PSF position ($\mu$m)')
    cbar = fig5.colorbar(CS)
    cbar.ax.set_ylabel('Normalized intensity')
if PSF==0:
    ax2.set_title('Contrast with z-position')
    ax2.set_xlabel('Angle(°)')
    ax2.set_ylabel('Diffusors position ($\mu$m)')
    cbar = fig5.colorbar(CS)
    cbar.ax.set_ylabel('Normalized intensity')
    
if choose_data==1:
    fig5.savefig('z_scan_{0}_split_lens_{1}_PSF_{2}.png'.format(name_file,split_lens,PSF), dpi=800, bbox_inches='tight', pad_inches=0) # For illustration            
if choose_data==0:
    fig5.savefig('z_scan_theo_h_{0}_w_{1}_PSF_{2}.png'.format(lens_heigth,lens_period,PSF), dpi=800, bbox_inches='tight', pad_inches=0) # For illustration           