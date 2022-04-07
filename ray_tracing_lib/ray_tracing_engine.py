from ray_tracing_lib.libraries_import import *
from ray_tracing_lib.general_functions import *

def Ray_tracing_engine(F,theta_lim,inc_x,nb_ray_calc_max,color_def,X_lens,Y_lens,step_calc,n_PC,n_air,y_lim_top):
    
    n=n_air/n_PC
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
    return out_rays