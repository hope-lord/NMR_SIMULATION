# packages
from math import pi,cos, sin
from RK4 import rk4
from FFT import FFT

#differential equation governing the amplitude of spin 
def diff_eq_of_spin(y,t,args):
    w0 = args[0]
    w = args[1]
    w1 = args[2]
    ar,ai,br,bi = y[0],y[1],y[2],y[3]
    dar = 0.5*(w0*ai+w1*cos(w*t)*bi)
    dai = -0.5*(w0*ar+w1*cos(w*t)*br)
    dbr = -0.5*(w0*bi-w1*cos(w*t)*ai)
    dbi = 0.5*(w0*br-w1*cos(w*t)*ar)
    return [dar,dai,dbr,dbi]

# probability that the spin is downward
def prob_down(dt,N,w,w1,w0=1e9):
    t, y = rk4(diff_eq_of_spin,[1,0,0,0],0,dt*N,N,args=[w0,w,w1])
    pb = []
    for i,j in zip(y[2],y[3]):
        pb.append(i**2+j**2)
    return t,pb
# power spectrum obtained by fourier transform of flipping probability
def power_spectrum(dt,N,w,w1,w0=1e9):
    t,pb = prob_down(dt,N,w,w1,w0)
    yw = FFT(pb)
    Yw = []
    wp = []
    for i in range(N):
        wp.append(i*2*pi/dt/(N-1))
        Yw.append(abs(yw[i])**2)
    return wp , Yw
 

