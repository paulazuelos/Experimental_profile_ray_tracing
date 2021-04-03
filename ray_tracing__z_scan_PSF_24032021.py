import math
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from scipy.interpolate import griddata
from scipy.interpolate import splev, splrep
from scipy import signal
from scipy.spatial import distance
import os

# P.Azuelos Idemia R&T team 2020
# Forward ray tracing algorithm;  LQA project 
# study of the effect of the lens shape
# study of the effect of the engraving and color matrix position
# render above lenses with an angle
# compare sperical and parabolic lens profile

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
def theoretical_lens_profile(a0=90,h0=12,n=1.5,nb_lens=20,split_lens=45):
    # Theoretical lens profil to test the program
    
    
    X_profile=[]
    Y_profile=[]
    for i in range(0,nb_lens):
    # geometrical parameters for lens i      
    
        a=a0 # diameter of the lens base
        h=h0 # Lens height 
        R = ((a/2)**2+h**2)/(2*h) # Curvature radius
        F=n*R/(n-1) # Focal length
        X_0 = 0 # Lens center
        Y_0 = 0-R+h # Lens center
        
        Theta_1=math.acos(((a/2)-X_0)/R)*180/math.pi
        Theta_2=math.acos((-(a/2)-X_0)/R)*180/math.pi
    
        Theta_list =np.linspace(Theta_1, Theta_2, num=split_lens) # Angles in the lens curvature
        # Lens profile
        actual_lens_X= R* np.cos (Theta_list* math.pi/ 180) + X_0+a*i+0.001*i
        actual_lens_X=np.flip(actual_lens_X)
        X_profile=np.concatenate((X_profile, actual_lens_X))
        actual_lens_Y= R* np.sin (Theta_list* math.pi/ 180) + Y_0
        actual_lens_Y=np.flip(actual_lens_Y)
        Y_profile=np.concatenate((Y_profile, actual_lens_Y))
    return(X_profile, Y_profile,F)
    
def theta_vect(x1, y1, x2, y2):
    len1 = (x1*x1+y1*y1)**(1/2)
    len2 = (x2*x2+y2*y2)**(1/2)
    x1 = x1 / len1
    y1 = y1 / len1
    x2 = x2 / len2
    y2 = y2 / len2

    theta = (np.arctan2(y2, x2) - np.arctan2(y1, x1)) * 180 / np.pi
    return theta

def snell_descarte(n0, n1, theta0):
    theta1=np.arcsin((n0 / n1) * np.sin(theta0 * np.pi / 180)) * (180 / np.pi)
    return theta1
#test =theta_vect(1,1,0,1)
def n_PC_func(wavelength):

    return (1+1.4182/(1-0.021304/wavelength**2))**.5  
# https://refractiveindex.info/?shelf=organic&book=polycarbonate&page=Sultanova
   
#%% input 
    
####### input parameters #####################################################

   
n_PC=n_PC_func(0.6) # refractive index of the polycarbonate

n_air=1 # refractive index of the air
n=n_air/n_PC

### input variable for the theoretical lens profile (choose_data=0)

lens_period=90 # pitch between two lenses
lens_heigth=17 # heigth of the lens (cylindrical profile)
nb_lens=10 # number of lens to be computed

### input variable for the experimental lens profile (choose_data=1)

expected_width=140 # used for smoothing and peak finding of the raw experimental datas
split_lens=math.floor(expected_width/4) # discretization steps per lens (typically between 90 and 200) check the results for speed optimisation
filter_detect_peak=0.2

directory="C:/Users/paulazue/OneDrive - myidemia/Doc Idemia Azuelos Paul/R&D/Lens related subjects/CO2_microlens/Workshop_CO2_nov_2020/Ray tracing simulation"
os.chdir(directory)
VTT='SLI profile 2.txt' 
FourPlate='125005 3.8.txt'
CO2='CO2 lens profiles/'
name_file='C39P3 HRX_MM  - Extracted profile.txt'

x_crop_min=100
x_crop_max=15000
# variables to drive the software parts

choose_data=1# 0 generate a theoretical lens array, 1 interpolate an experimental lens array, other number don't calculate again
lamination_plate=0 # sometimes we want to check the difference between plates profile and actual laminated lens (it just flip the data along y-axis)
calculation=1 # calculate only when the lens profile and the calculation points are set up
plotting_experimental=1
ray_tracing=1
# Angle of calculation
theta_lim=2*math.asin(n_air/n_PC)*180/math.pi
out_angle_lim=90 # range of angles for the plot
split_out_angle=2*out_angle_lim # for the discretization of angle range
range_angle=np.linspace(-out_angle_lim,out_angle_lim,split_out_angle)
# PSF calculation
PSF=1
nb_step=25
z_depth=np.linspace(150,500,nb_step)
Z_scan_PSF=np.zeros((range_angle.size,nb_step))



#%% Calculation and interpolation of the lens profile

###### theoretical spherical lens profile ######
if choose_data==0:
    expected_width=lens_period
    [X_lens0, Y_lens_ini,F]=theoretical_lens_profile(lens_period,lens_heigth,n_PC,nb_lens,split_lens)
    x_grid=np.linspace(np.min(X_lens0),np.max(X_lens0),X_lens0.size)
    X_lens=x_grid
    X_lens0=np.array(X_lens0)
    Y_lens_ini=np.array(Y_lens_ini)
    
    tck = splrep(X_lens0, Y_lens_ini)
    Y_lens0 = splev(X_lens, tck)
    
    dx=X_lens0[1]-X_lens0[0] # experimental incrementation along x-direction (expected fixed)
    
    target_width=round(expected_width/dx) # estimation of the number of measurements per period
    # Apply a moving average filter to find the peaks
    x_values = np.linspace(-1, 1, math.floor(target_width))
    weight=gaussian(x_values, 0, 0.04)
    average_filter=weight/np.sum(weight)
    Ylens_peak=np.correlate(Y_lens0,average_filter, "same")
    
    peakind = signal.find_peaks_cwt(Ylens_peak, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
    # remove the first and last peak in order to keep only full period
    peakind= np.delete(peakind, 0, 0)
    peakind= np.delete(peakind, -1, 0)
    # 
    lens_period=np.zeros((peakind.size,1))
    for i in range(1,peakind.size):
        lens_period[i-1]=(peakind[i]-peakind[i-1])*dx
    lens_period=lens_period[~np.all(lens_period == 0, axis=1)]
    lens_period_mean=np.mean(lens_period)
    
    ############ show detection of lens position ############
    if plotting_experimental==1:
        
        f1 = plt.figure()
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
        plt.plot(X_lens0, Y_lens0, '.')
        plt.plot(X_lens0, Ylens_peak, '.')
        plt.plot(X_lens0[peakind], Ylens_peak[peakind], '*')
        
        
        x1,x2,y1,y2 = plt.axis()
        #plt.axis((T_grid_min,T_grid_max,y1,90))
        plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
        
        plt.show()
    
    ######################################################
    
    
    hollow_ind = signal.find_peaks_cwt(-Y_lens0, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
    
    hollow_ind= np.delete(hollow_ind, 0, 0)
    hollow_ind= np.delete(hollow_ind, -1, 0)
    nb_lens_real=len(hollow_ind)
    print(hollow_ind)
    
    ############ show detection of hollow position ############
    
    if plotting_experimental==1:
        
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
        plt.plot(X_lens0, Y_lens0, '.')
        #plt.plot(X_lens0, Ylens_peak, '.')
        plt.plot(X_lens0[hollow_ind], Y_lens0[hollow_ind], '*')
        
        
        x1,x2,y1,y2 = plt.axis()
        #plt.axis((T_grid_min,T_grid_max,y1,90))
        plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
        
        plt.show()
    
    #####################################################################
    ######################### Apply the interpolation ############################
   
    x_grid=X_lens0

    Y_lens = Y_lens0
    Y_lens_grad_1=np.gradient(Y_lens)
    Y_lens_grad_2=np.gradient(Y_lens_grad_1)
    dx=X_lens[1]-X_lens[0]
    target_width=split_lens
    lens_period_mean=split_lens*dx
    # Apply a filter to detect the peaks
    x_values = np.linspace(-1, 1, math.floor(target_width))
    weight=gaussian(x_values, 0, 0.2)
    average_filter=weight/np.sum(weight)
    Ylens_peak=np.correlate(Y_lens,average_filter, "same")
    
    peakind = signal.find_peaks_cwt(Ylens_peak, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
    
    
    if plotting_experimental==1:
        
       
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
        plt.plot(X_lens, Ylens_peak, '.')
        #plt.plot(X_lens0, Ylens_peak, '.')
        plt.plot(X_lens[peakind], Ylens_peak[peakind], '*')
        
        x1,x2,y1,y2 = plt.axis()
        #plt.axis((T_grid_min,T_grid_max,y1,90))
        plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
        
        plt.show()
    
    hollow_ind = signal.find_peaks_cwt(-Ylens_peak, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
    # estimate the height of the lens
    mean_heigth=np.mean(abs(Y_lens[peakind[0:min(peakind.size,hollow_ind.size)]]-Y_lens[hollow_ind[0:min(peakind.size,hollow_ind.size)]]))
    # recentered the datas at y=0
    Y_lens =Y_lens-np.amin(Y_lens)
    # estimate the filling factor
    above_height=np.count_nonzero(Y_lens > mean_heigth/2)
    below_height=len(Y_lens)-above_height
    FF=abs(below_height-above_height)/(len(Y_lens))+0.5 ## estimate the filling factor of the lens
    # estimate the lens radius
    R = (((lens_period_mean*FF)/2)**2+mean_heigth**2)/(2*mean_heigth) # Curvature radius
    F=n_PC*R/(n_PC-1) # Focal length
    # estimate the lens pitch
    lens_period=round(lens_period_mean)
###### experimental lens profile ######
if choose_data==1:

    data_input=np.char.replace(np.genfromtxt(CO2+name_file,dtype='str'), ',', '.').astype(float)
    data_0=data_input[x_crop_min:x_crop_max,:]
    X_lens0=data_0[:,0]
    if lamination_plate==1:
        Y_lens0=-data_0[:,1]+np.amin(-data_0[:,1])
    if lamination_plate==0:
        Y_lens0=data_0[:,1]
    
    dx=X_lens0[1]-X_lens0[0] # experimental incrementation along x-direction (expected fixed)
    
    target_width=round(expected_width/dx) # estimation of the number of measurements per period
    # Apply a moving average filter to find the peaks
    x_values = np.linspace(-1, 1, math.floor(target_width))
    weight=gaussian(x_values, 0, filter_detect_peak)
    average_filter=weight/np.sum(weight)
    Ylens_peak=np.correlate(Y_lens0,average_filter, "same")
    
    peakind = signal.find_peaks_cwt(Ylens_peak, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
    # remove the first and last peak in order to keep only full period
    peakind= np.delete(peakind, 0, 0)
    peakind= np.delete(peakind, -1, 0)
    # 
    lens_period=np.zeros((peakind.size,1))
    for i in range(1,peakind.size):
        lens_period[i-1]=(peakind[i]-peakind[i-1])*dx
    lens_period=lens_period[~np.all(lens_period == 0, axis=1)]
    lens_period_mean=np.mean(lens_period)
    
    ############ show detection of lens position ############
    if plotting_experimental==1:
        
        f1 = plt.figure()
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
        plt.plot(X_lens0, Y_lens0, '.')
        plt.plot(X_lens0, Ylens_peak, '.')
        plt.plot(X_lens0[peakind], Ylens_peak[peakind], '*')
        
        
        x1,x2,y1,y2 = plt.axis()
        #plt.axis((T_grid_min,T_grid_max,y1,90))
        plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
        
        plt.show()
    
    ######################################################
    
    
    hollow_ind = signal.find_peaks_cwt(-Y_lens0, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
    
    hollow_ind= np.delete(hollow_ind, 0, 0)
    hollow_ind= np.delete(hollow_ind, -1, 0)
    nb_lens_real=len(hollow_ind)
    
    
    ############ show detection of hollow position ############
    
    if plotting_experimental==1:
        
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
        plt.plot(X_lens0, Y_lens0, '.')
        #plt.plot(X_lens0, Ylens_peak, '.')
        plt.plot(X_lens0[hollow_ind], Y_lens0[hollow_ind], '*')
        
        
        x1,x2,y1,y2 = plt.axis()
        #plt.axis((T_grid_min,T_grid_max,y1,90))
        plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
        
        plt.show()
    
    #####################################################################
    ######################### Apply the interpolation ############################
    nb_point=(nb_lens_real-1)*split_lens
    x_grid=np.linspace(X_lens0[hollow_ind[0]],X_lens0[hollow_ind[nb_lens_real-1]],nb_point)
    X_lens=x_grid
    X_lens0=np.array(X_lens0)
    Y_lens0=np.array(Y_lens0)
    tck = splrep(X_lens0, Y_lens0)
    Y_lens = splev(X_lens, tck)
    Y_lens_grad_1=np.gradient(Y_lens)
    Y_lens_grad_2=np.gradient(Y_lens_grad_1)
    dx=X_lens[1]-X_lens[0]
    target_width=split_lens
    lens_period_mean=split_lens*dx
    # Apply a filter to detect the peaks
    x_values = np.linspace(-1, 1, math.floor(target_width))
    weight=gaussian(x_values, 0, 0.2)
    average_filter=weight/np.sum(weight)
    Ylens_peak=np.correlate(Y_lens,average_filter, "same")
    
    peakind = signal.find_peaks_cwt(Ylens_peak, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
    
    
    
    if plotting_experimental==1:
        
       
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
        plt.plot(X_lens, Ylens_peak, '.')
        #plt.plot(X_lens0, Ylens_peak, '.')
        plt.plot(X_lens[peakind], Ylens_peak[peakind], '*')
        
        x1,x2,y1,y2 = plt.axis()
        #plt.axis((T_grid_min,T_grid_max,y1,90))
        plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
        
        plt.show()
    
    hollow_ind = signal.find_peaks_cwt(-Ylens_peak, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
    # estimate the height of the lens
    mean_heigth=np.mean(abs(Y_lens[peakind[0:min(peakind.size,hollow_ind.size)]]-Y_lens[hollow_ind[0:min(peakind.size,hollow_ind.size)]]))
    # recentered the datas at y=0
    Y_lens =Y_lens-np.amin(Y_lens)
    # estimate the filling factor
    above_height=np.count_nonzero(Y_lens > mean_heigth/2)
    below_height=len(Y_lens)-above_height
    FF=abs(below_height-above_height)/(len(Y_lens))+0.5 ## estimate the filling factor of the lens
    # estimate the lens radius
    R = (((lens_period_mean*FF)/2)**2+mean_heigth**2)/(2*mean_heigth) # Curvature radius
    F=n_PC*R/(n_PC-1) # Focal length
    # estimate the lens pitch
    lens_period=round(lens_period_mean)
    
for idx, nb_step_calculation_points in enumerate(z_depth):
    ### input variables for the positionning and discretization of the calculation points

    if ray_tracing==1:
        step_calc=5
    if ray_tracing==0:
        step_calc=1
    y_lim_calcualtion=-nb_step_calculation_points
    y_lim_top=math.floor(nb_step_calculation_points/4)
    nb_points_max=100000
    nb_cut=8
    
    ### input variables for calculation
    
    nb_ray_calc_max=800000
    
    ### input variables for plot
    
    colors = np.array([(1, 1, 0), (0, 1, 1), (1, 0, 1),(0.5,0.5,0.5)])
    
    #############################################################################
    #############################################################################
        
    #%% Definition of the diffusion points (matrix or engraving position)
    # the points can be positionned anywhere below the lens 
    
    y_grid=np.linspace(0,y_lim_calcualtion,nb_step_calculation_points+1)
    
    inc_x=abs(x_grid[0]-x_grid[1])
    inc_y=abs(y_grid[0]-y_grid[1])
    
    pitch_matrix=math.floor(2*lens_period/nb_cut)
    color_def=np.zeros((nb_points_max,5))
    count_color=0
    # here I parametrize a Lasink matrix with 4 colors
    # when we study the engraving, we can choose to work only with one color 
    
    spot_size=1
    mini=np.ones((peakind.size,1))*50
    for i in range(0,x_grid.size-1):
        count_color=count_color+1
        for count in range(0,peakind.size):
            if abs(x_grid[i]-x_grid[peakind[count]])<spot_size and abs(x_grid[i]-x_grid[peakind[count]])<mini[count]:
                mini[count]=abs(x_grid[i]-x_grid[peakind[count]])
                

    if PSF==1:
        
        count_color=0          
        
        for i in range(0,x_grid.size-1):
            count_color=count_color+1
            for count in range(0,peakind.size):
                if abs(x_grid[i]-x_grid[peakind[count]])==mini[count]:
                    color_def[count_color,0:3]=[2, 0, 0]
                    color_def[count_color,3]=x_grid[i]
                    color_def[count_color,4]=y_grid[-1]
                
    if PSF==0:   
        
        step_central_color=25
        step=math.floor(step_central_color/inc_x)
        """
        for i in range(0,x_grid.size-1):
            count_color=count_color+1
            
            color_def[count_color,0:3]=[1, 0, 0]
            color_def[count_color,3]=x_grid[i]
            color_def[count_color,4]=y_grid[-1]
        """    
        count_color=0  
         
        for i in range(0,x_grid.size-1):
            center=0
            for count in range(0,peakind.size):
                count_color=count_color+1
                
                if abs(x_grid[i]-x_grid[peakind[count]])< step_central_color:
                    
                    color_def[count_color,0:3]=[2, 0, 0]
                    color_def[count_color,3]=x_grid[i]
                    color_def[count_color,4]=y_grid[-1]
                    center=1
            if center==0:
                
                color_def[count_color,0:3]=[1, 0, 0]
                color_def[count_color,3]=x_grid[i]
                color_def[count_color,4]=y_grid[-1]
                
    
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
    """
    for peak in peakind:
        plt.plot(x_vect*X_lens[peak],y_grid,'-r')
    """    
    plt.scatter(color_def[:,3], color_def[:,4], s=1, c=colors[color_matrix])
    
    x1,x2,y1,y2 = plt.axis()
    #plt.axis((T_grid_min,T_grid_max,y1,90))
    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
    f2.savefig('measured_conf_{0}.png'.format(name_file), dpi=800, bbox_inches='tight', pad_inches=0)
    plt.show()
    
    #%% Calculation (forward ray tracing vectors calculation)
    
    if calculation==1:
        ##### Calculate all the forward rays angles direction
        
        n_inc_calc=math.floor((F*math.tan((theta_lim)*math.pi/180))/inc_x)   # Angles to consider
        #n_inc_calc=40*inc_x
        out_rays=np.zeros((nb_ray_calc_max,12)) # pre-allocated memory the size can be changed if memory issue appears
        inc_color=np.ones((1,12))
        values_test=np.ones((1,11))
        inc_test=np.ones((1,11))
        count=0
        actual_y=-10000
        
        for j in range(0,color_def[:,1].size):
            if color_def[j,0]==1:
                n_inc_calc=math.floor((F*math.tan((40)*math.pi/180))/inc_x)
                test=np.where((X_lens[:]-color_def[j,3])**2<n_inc_calc**2,0,1) # select only the the lens profile portion around the point to be ray traced
            if color_def[j,0]==2:
                n_inc_calc=math.floor((F*math.tan((40)*math.pi/180))/inc_x)
                test=np.where((X_lens[:]-color_def[j,3])**2<n_inc_calc**2,0,1) # select only the the lens profile portion around the point to be ray traced
            select_x=X_lens[test==0]
            select_y=Y_lens[test==0]
            for i in range(0,select_y.size-1,step_calc):
                
                #if abs(actual_y-select_y[i])>inc_y:
                actual_y=select_y[i]
                # 0 deg vector
                x_0=0
                y_0=1
                # color matrix to lens surface vector
                x_diff_surf=color_def[j,3]-select_x[i]
                y_diff_surf=select_y[i]-color_def[j,4]
                # color matrix to lens surface vector
                x_diff_surf_1=color_def[j,3]-select_x[i+1]
                y_diff_surf_1=select_y[i+1]-color_def[j,4]
                
                # perpendicular to the surface vector
                x_diff_perp_surf=select_y[i-1]-select_y[i]
                y_diff_perp_surf=select_x[i]-select_x[i-1]
                
                theta0=theta_vect(x_0,y_0,x_diff_surf,y_diff_surf) # angle between y-axis and diffused ray
                theta1=theta_vect(x_0,y_0,x_diff_surf_1,y_diff_surf_1)
                theta_perp=-theta_vect(x_0,y_0,x_diff_perp_surf,y_diff_perp_surf) # angle between y-axis and vector perpendicular to the lens surface
                theta_surf=-(theta0-theta_perp) # angle between perpendicular to the lens surface and input diffused ray
            
                if abs(theta_surf)<theta_lim:
                    
                    theta_out=snell_descarte(n_PC,n_air,theta_surf) # final output angle between the output ray and the y axis
                    theta_fin=theta_perp-theta_out
                    # Calculations of the reflectivity and the transmission at the lens interface
                    Rperpendicular=((n*math.cos(theta_surf*math.pi/180)-math.cos(theta_out*math.pi/180))/(n*math.cos(theta_surf*math.pi/180)+math.cos(theta_out*math.pi/180)))**2
                    Rparallel=((math.cos(theta_surf*math.pi/180)-n*math.cos(theta_out*math.pi/180))/(math.cos(theta_surf*math.pi/180)+n*math.cos(theta_out*math.pi/180)))**2
                    T=1-0.5*(Rperpendicular+Rparallel)
                    elem_surface=math.sqrt(abs(select_y[i]-select_y[i+1])**2+abs(select_x[i]-select_x[i+1])**2)
                    x_top_ray=select_x[i]+abs(select_y[i]-y_lim_top)*math.tan(theta_fin*math.pi/180)
                    delta_theta=abs(theta0-theta1)
                    # registration only if the rays are not internally reflected
                    count=count+1
                    inc_color[0,0]=theta0
                    inc_color[0,1]=theta_surf
                    inc_color[0,2]=theta_fin
                    inc_color[0,3]=color_def[j,3]
                    inc_color[0,4]=color_def[j,0]
                    inc_color[0,5]=select_x[i]
                    inc_color[0,6]=delta_theta
                    inc_color[0,7]=math.cos(theta0*math.pi/180)*T*math.cos(theta_surf*math.pi/180)*delta_theta # the formulae of the intensity may change in the future 
                    inc_color[0,8]=x_top_ray
                    inc_color[0,9]=color_def[j,4]
                    inc_color[0,10]=select_y[i]
                    inc_color[0,11]=y_lim_top
                    
                    out_rays[count,:]=inc_color
        out_rays=out_rays[~np.all(out_rays == 0, axis=1)]
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
        if choose_data==1: 
            """
            color1=np.sum(value_mean[1:nb_lens_real-2,:,0],axis=0)
            color2=np.sum(value_mean[1:nb_lens_real-2,:,1],axis=0)
            
            contraste=(color2-color1)/(color2+color1)
            contraste = np.nan_to_num(contraste) 
            """
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
            
            if PSF==0:
                
                mean_2=np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)
                mean_1=np.sum(value_mean[1:nb_lens_real-2,:,1],axis=0)
                contrast=(mean_2-mean_1)/(mean_2+mean_1)
                
                Z_scan_PSF[:,idx]=contrast
            
            if PSF==1:
                Z_scan_PSF[:,idx]=np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)/max(np.sum(np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)),1)
            
            f1.savefig('x_theta_render_measured_{0}_split_lens_{1}.png'.format(name_file,split_lens).format(name_file), dpi=800, bbox_inches='tight', pad_inches=0) # For illustration
            color_matrix=np.array(color_def[:,0].astype(int))
            
            """
            f2 = plt.figure()
            Font_SIZE = 15
            plt.rc('font', size=Font_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=Font_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=Font_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=Font_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=Font_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=Font_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=Font_SIZE)  # fontsize of the figure title
            
            plt.xlabel('Angle (°)')
            plt.ylabel('Contraste between two colors')
    
            plt.plot(range_angle, contraste, '-r')
            f2.savefig('contraste.png', dpi=800, bbox_inches='tight', pad_inches=1) # For illustration
            color_matrix=np.array(color_def[:,0].astype(int))
            """
            X=X_lens[hollow_ind[1]:hollow_ind[2], np.newaxis]
            Y=Y_lens[hollow_ind[1]:hollow_ind[2], np.newaxis]
            
            A = np.hstack([X**2, X * Y, Y**2, X, Y])
            b = np.ones_like(X)
            x = np.linalg.lstsq(A, b)[0].squeeze()
            
            x_coord = np.linspace(X_lens[hollow_ind[1]],X_lens[hollow_ind[2]],100)
            y_coord = np.linspace(np.amin(Y_lens[hollow_ind[1]:hollow_ind[2]]),np.amax(Y_lens[hollow_ind[1]:hollow_ind[2]]),100)
            X_coord, Y_coord = np.meshgrid(x_coord, y_coord)
            Z_coord = x[0] * X_coord ** 2 + x[1] * X_coord * Y_coord + x[2] * Y_coord**2 + x[3] * X_coord + x[4] * Y_coord
            
    
            f3 = plt.figure()
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
            plt.plot(X_lens[hollow_ind[1]:hollow_ind[2]], Y_lens[hollow_ind[1]:hollow_ind[2]], '.-', label='measured_datas')
            CS=plt.contour(X_coord, Y_coord, Z_coord, levels=[1], colors=('r'), linestyle='--', linewidths=0)
            data_fit= CS.allsegs[0][0]
            data_exp=np.zeros((X_lens[hollow_ind[1]:hollow_ind[2]].shape[0],2))
            data_exp[:,0]=X_lens[hollow_ind[1]:hollow_ind[2]]
            data_exp[:,1]=Y_lens[hollow_ind[1]:hollow_ind[2]]
            distance_found=distance.cdist(data_exp,data_fit).min(axis=1)
            plt.plot(data_fit[:,0],data_fit[:,1], '--', label='ellipsoïd fit')
            plt.axis()
            plt.legend()
            #plt.axis((T_grid_min,T_grid_max,y1,90))
            plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
            f3.savefig('lens_profile_{0}_split_lens_{1}.png'.format(name_file,split_lens), dpi=800, bbox_inches='tight', pad_inches=0)
            plt.show()
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
                theta_one_lens=0.8*(math.pi/2-math.atan(abs(y_lim_calcualtion)/(lens_period_mean/2)))*(180/math.pi)
                for i in range(0,out_rays.shape[0]):
                    if out_rays[i,8]<X_lens[-1] and out_rays[i,8]>X_lens[0] and abs(out_rays[i,0])<theta_one_lens:
                        
                        plt.plot([out_rays[i,3],out_rays[i,5]], [out_rays[i,9],out_rays[i,10]],'-b', linewidth=0.25)
                        plt.plot([out_rays[i,5],out_rays[i,8]], [out_rays[i,10],out_rays[i,11]],'-b', linewidth=0.25)
                    
                x1,x2,y1,y2 = plt.axis()
                #plt.axis((T_grid_min,T_grid_max,y1,90))
                plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
                f4.savefig('measured_conf_ray_tracing_{0}_split_lens_{1}.png'.format(name_file,split_lens), dpi=800, bbox_inches='tight', pad_inches=0)
                plt.show()
            
            
        if choose_data==0:
            
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
            
            if PSF==0:
                
                mean_2=np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)
                mean_1=np.sum(value_mean[1:nb_lens_real-2,:,1],axis=0)
                contrast=(mean_2-mean_1)/(mean_2+mean_1)
                
                Z_scan_PSF[:,idx]=contrast
            
            if PSF==1:
                Z_scan_PSF[:,idx]=np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)/max(np.sum(np.sum(value_mean[1:nb_lens_real-2,:,2],axis=0)),1)
            
            
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
                theta_one_lens=0.8*(math.pi/2-math.atan(abs(y_lim_calcualtion)/(lens_period_mean/2)))*(180/math.pi)
                for i in range(0,out_rays.shape[0]):
                    if out_rays[i,8]<X_lens[-1] and out_rays[i,8]>X_lens[0] and abs(out_rays[i,0])<theta_one_lens:
                        
                        plt.plot([out_rays[i,3],out_rays[i,5]], [out_rays[i,9],out_rays[i,10]],'-b', linewidth=0.25)
                        plt.plot([out_rays[i,5],out_rays[i,8]], [out_rays[i,10],out_rays[i,11]],'-b', linewidth=0.25)
                    
                x1,x2,y1,y2 = plt.axis()
                #plt.axis((T_grid_min,T_grid_max,y1,90))
                plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
                f4.savefig('measured_conf_ray_tracing_spherical_profil_split_lens_{0}.png'.format(split_lens), dpi=800, bbox_inches='tight', pad_inches=0)
                plt.show()
            """
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
            plt.scatter(color_def[:,3], color_def[:,4], s=1, c=colors[color_matrix])
            
            color1=np.sum(value_mean[1:nb_lens-2,:,1],axis=0)
            color2=np.sum(value_mean[1:nb_lens-2,:,0],axis=0)
            
            contraste_normalized=(color2-color1)/(color2+color1)
            contraste=(color2-color1)
            contraste = np.nan_to_num(contraste) 
            
            
            f3 = plt.figure()
            Font_SIZE = 15
            plt.rc('font', size=Font_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=Font_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=Font_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=Font_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=Font_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=Font_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=Font_SIZE)  # fontsize of the figure title
            
            plt.xlabel('Angle (°)')
            plt.ylabel('Contraste between two colors')
            
            plt.plot(range_angle, contraste, '-r')
            f2.savefig('contraste.png', dpi=800, bbox_inches='tight', pad_inches=1) # For illustration
            color_matrix=np.array(color_def[:,0].astype(int))
            int_col1=np.sum(np.sum(value_mean[1:nb_lens-2,28:46,1],axis=0),axis=0)
            int_col2=np.sum(np.sum(value_mean[1:nb_lens-2,28:46,0],axis=0),axis=0)
            """
#%%            
z_depth_map,range_angle_map=np.meshgrid(z_depth,range_angle)

#%%

fig5, ax2 = plt.subplots(constrained_layout=True)
origin = 'lower'

CS = ax2.contourf(range_angle_map, z_depth_map,Z_scan_PSF, 100, cmap=plt.cm.bwr, origin=origin)

if PSF==1:
    focal_pos_z=z_depth_map[np.unravel_index(np.argmax(Z_scan_PSF, axis=None), Z_scan_PSF.shape)]
    focal_pos_theta=range_angle_map[np.unravel_index(np.argmax(Z_scan_PSF, axis=None), Z_scan_PSF.shape)]
    ax2.scatter(focal_pos_theta, focal_pos_z,s=500,marker="x", color="black")
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