from libraries_import import *
"""
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
    f2.savefig('measured_conf_{0}.png'.format(name_file), dpi=800, bbox_inches='tight', pad_inches=0)
    plt.show()
"""