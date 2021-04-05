from libraries_import import *

def Emitters_positionning(y_lim_calculation,nb_step_calculation_points,x_grid,y_grid,lens_period,nb_cut,peakind,nb_points_max,PSF,dist_center):

    y_grid=np.linspace(0,y_lim_calculation,nb_step_calculation_points+1)
    
    inc_x=abs(x_grid[0]-x_grid[1])
    inc_y=abs(y_grid[0]-y_grid[1])
    dist_center_norm=math.floor(dist_center/inc_x)*inc_x
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
                    color_def[count_color,3]=x_grid[i]+dist_center_norm
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
    return color_def