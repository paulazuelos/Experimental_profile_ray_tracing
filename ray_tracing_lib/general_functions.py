from ray_tracing_lib.libraries_import import *

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
   