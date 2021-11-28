from pylab import*
from matplotlib.widgets import Slider, Button
from nmr import rk4
#from scipy import integrate as inte

def slider(l,initcond,pulse,w0p,wp,w1p):
    axomega = plt.axes([0.25, 0.15, 0.65, 0.03])
    axomega1 = plt.axes([0.25, 0.1, 0.65, 0.03])
    axomega0 = plt.axes([0.25, 0.05, 0.65, 0.03])

    omega = Slider(axomega, '$\omega$', 20*pi, 130*pi, w0p,valstep=pi)
    omega0 = Slider(axomega0, '$\omega_0$', 20*pi, 130*pi, wp,valstep=pi)
    omega1 = Slider(axomega1, '$\omega_1$', 0.0,130*pi, w1p, valstep=pi)


    def update(val):
        w = omega.val
        w1 = omega1.val
        w0 = omega0.val
        t,M=rk4(dM_dt_rot,initcond,0,2,3000,[w0,w,w1,1,0.3,pulse])
    
        l.set_data_3d(*M)
 
    # Call update function when slider value is changed
    omega.on_changed(update)
    omega1.on_changed(update)
    omega0.on_changed(update)



##############################################################
##############################################################

def dM_dt_rot(M,t,args): # for rotating frame
    
    M0 = 1
    w0 = args[0]
    w = args[1]
    w1 = args[2]
    T1 = args[3]
    T2 = args[4]
    delta = w-w0
    tu = 0
    if args[5]!=0:
        pulse_type = args[5]*pi/180
        tu =  pulse_type/(delta**2+w1**2)**0.5
     
    w1 = args[2] if t < tu else 0
    dM_x = -M[0]/T2 + delta*M[1]#  -M[2]*w1*cos(w*t)
    dM_y = -delta*M[0]-M[1]/T2 + M[2]*w1#*sin(w*t)
    dM_z = +(M0-M[2])/T1 -w1*M[1]#*cos(w*t) #+w1*M[0]*sin(w*t)
    return [dM_x, dM_y,dM_z]
#exit()



################################
################################
################################
