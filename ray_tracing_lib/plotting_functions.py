from ray_tracing_lib.libraries_import import *

def plot_profile(X,Y):
    
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
    plt.plot(X, Y, '.')

    x1,x2,y1,y2 = plt.axis()
    #plt.axis((T_grid_min,T_grid_max,y1,90))
    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
    
    plt.show()

