from ray_tracing_lib.libraries_import import *
from ray_tracing_lib.general_functions import *


def Interp_profile_theo(lens_period,lens_heigth,n_PC,nb_lens,split_lens):
    
    """
    
    This function interpolate a theoretical lens profile which can be used as a test profile 
    for the positionning of the emitters and the ray tracing engine
    
    Args:
        lens_period : period of the lens array
        lens_heigth : height of the cylindrical lens
        n_PC : refractive index of the polycarbonate material
        nb_lens : number of lenses
        split_lens : 
            
    Returns:
        X_lens0: X_coordinates of the lens profile
        Y_lens0: Y_coordinates of the lens profile
        F: Focal length of the lens profile
        Y_lens: interpolated Y_coordinates of the lens profile
        peakind: x coordinates of the maximas

    """
    expected_width=lens_period
    [X_lens0, Y_lens0,F]=theoretical_lens_profile(lens_period,lens_heigth,n_PC,nb_lens,split_lens)
    x_grid=np.linspace(np.min(X_lens0),np.max(X_lens0),X_lens0.size)
    X_lens=x_grid
    X_lens0=np.array(X_lens0)
    Y_lens0=np.array(Y_lens0)
    
    tck = splrep(X_lens0, Y_lens0)
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
    
    hollow_ind = signal.find_peaks_cwt(-Y_lens0, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
    
    hollow_ind= np.delete(hollow_ind, 0, 0)
    hollow_ind= np.delete(hollow_ind, -1, 0)
    nb_lens_real=len(hollow_ind)

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
    return X_lens, Y_lens0,F,Y_lens,peakind,nb_lens_real,lens_period_mean

def Interp_profile_exp(file_path,x_crop_min,x_crop_max,expected_width,filter_detect_peak,lamination_plate,split_lens,n_PC):
        
        """
        
        This function interpolate an experimental lens profile which can be used
        for the positionning of the emitters and the ray tracing engine
        
        Args:
            lens_period : period of the lens array
            lens_heigth : height of the cylindrical lens
            n_PC : refractive index of the polycarbonate material
            nb_lens : number of lenses
            split_lens : 
                
        Returns:
            X_lens0: X_coordinates of the lens profile
            Y_lens0: Y_coordinates of the lens profile
            F: Focal length of the lens profile
            Y_lens: interpolated Y_coordinates of the lens profile
            peakind: x coordinates of the maximas
    
        """
        
        data_input=np.char.replace(np.genfromtxt(file_path,dtype='str'), ',', '.').astype(float)
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

        hollow_ind = signal.find_peaks_cwt(-Y_lens0, np.arange(target_width-target_width/1.25,target_width),noise_perc=20)
        
        hollow_ind= np.delete(hollow_ind, 0, 0)
        hollow_ind= np.delete(hollow_ind, -1, 0)
        nb_lens_real=len(hollow_ind)
        
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
        #plt.plot(X_lens0, Ylens_peak, '.')
        plt.plot(X_lens[hollow_ind], Y_lens[hollow_ind], '*')
        
        
        x1,x2,y1,y2 = plt.axis()
        #plt.axis((T_grid_min,T_grid_max,y1,90))
        plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
        
        plt.show()
    
    
        return X_lens, Y_lens0, F, Y_lens, peakind,nb_lens_real,lens_period_mean