from nmr import rk4


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




################################
################################
################################
