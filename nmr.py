# packages
from math import pi,cos, sin
## Runge-Kutta 4 algorithm 
def mul(a,lst): #multiplication of elements of a list with a scalar
    return [a*i for i in lst]
def add(*lst): #list addition similar to vector addition
    return [sum(i) for i in zip(*lst)]
def call_func(f,args,*x): #helps in passing the arguments
    if args==None:
        #print("Is not it")
        return f(*x)
    else: return f(*x,args)

def rk4(f,y0,x0,xstop,N,args=None):
    h=(xstop-x0)/(N-1)
    prev_result = y0
    y=[[prev_result[i]] for i in range(len(y0))]
    c = [x0]
    iter=0
    while iter<N-1:
        iter+=1
        #finding k1/2 and k1/6
        d = call_func(f, args, prev_result, x0)
        k1_by_2 = mul(h / 2, d)
        k1_by_6=mul(h / 6, d)
        #finding k2/2 and k2/3
        yn = add(prev_result, k1_by_2)
        d = call_func(f, args, yn, x0 + h / 2)
        k2_by_2=mul(h/2,d)
        k2_by_3=mul(h/3,d)
        #finding k3
        yn = add(prev_result, k2_by_2)
        d = call_func(f, args, yn, x0 + h / 2)
        k3 = mul(h, d)
        k3_by_3= mul(h/3, d)
        #finding k4
        yn = add(prev_result, k3)
        d = call_func(f, args, yn, x0+h)
        k4_by_6=mul(h/6,d)
        #finding the correct yn
        yn=add(prev_result, k1_by_6, k2_by_3, k3_by_3, k4_by_6)
        prev_result=yn
        for i in range(len(y0)): y[i].append(prev_result[i])
        x0 += h
        c.append(x0)
    return c,y

# Cooley-Tukey FAST FOURIER TRANSFORM
def FFT(x):
    N = len(x)
    
    if N == 1:
        return x
    elif N%2!=0:
        print("Length of the array must be 2^n")
        exit(1)
    else:
        x_even = []
        x_odd = []
        factor = []
        for i in range(N): #dividin into odd and even index
            if i%2==0: x_even.append(x[i]) 
            else : x_odd.append(x[i])
            factor.append( cos(2*pi*i/N)-1j*sin(2*pi*i/N)) #multiplicating factor
        X_even = FFT(x_even) # FFT of even part
        X_odd = FFT(x_odd) #FFT of odd part

        part1 = []
        part2 = []
        for i in range(N//2):
            part1.append(X_even[i]+factor[i]*X_odd[i])
            part2.append(X_even[i]+factor[N//2+i]*X_odd[i])
        return [*part1,*part2]

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


def prob_down(dt,N,w,w1,w0=1e9):
    t, y = rk4(diff_eq_of_spin,[1,0,0,0],0,dt*N,N,args=[w0,w,w1])
    pb = []
    for i,j in zip(y[2],y[3]):
        pb.append(i**2+j**2)
    return t,pb
def power_spectrum(dt,N,w,w1,w0=1e9):
    t,pb = prob_down(dt,N,w,w1,w0)
    yw = FFT(pb)
    Yw = []
    wp = []#linspace(0,1/dt,N)
    for i in range(N):
        wp.append(i*2*pi/dt/(N-1))
        Yw.append(abs(yw[i])**2)
    return wp , Yw
 

