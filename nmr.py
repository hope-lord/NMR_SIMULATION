from my_pkg.ode import rk4
import numpy as np
import matplotlib.pylab as plt
from scipy import integrate

def up_down_prob(y,t,args):
    w0 = args[0]
    w = args[1]
    w1 = args[2]
    ar,ai,br,bi = y[0],y[1],y[2],y[3]
    dar = 0.5*(w0*ai+w1*np.cos(w*t)*bi)
    dai = -0.5*(w0*ar+w1*np.cos(w*t)*br)
    dbr = -0.5*(w0*bi-w1*np.cos(w*t)*ai)
    dbi = 0.5*(w0*br-w1*np.cos(w*t)*ar)
    return [dar,dai,dbr,dbi]


######### UPDOWN PROBABILITY
t , y = rk4(up_down_prob,[1,0,0,0],0,2e-7,8000,args=[1e9,1e9,0.5e8])
t = np.array(t)
y = np.array(y)

#plt.plot(t*1e9,y[2]**2+y[3]**2,'g-')
#plt.xlim(0,1000)
#plt.ylim(0,1)
#plt.grid()
#plt.clf()
#exit()



######### FLIPPING PROBABILITY
'''
w0 = 1e9
w = np.linspace(0.9,1.1,1000)*1e9


def flip_prob(omegas,time_end,N):
    t , y = rk4(up_down_prob,[1,0,0,0],0,time_end,N,args=omegas)
    t = np.array(t)
    y = np.array(y)
    p_max = np.max(y[3]**2+y[2]**2)
    return p_max

pflip1 = []
w1 = 1e8
for i in w:  
    pmax = flip_prob([w0,i,w1],1e-7,500)
    pflip1.append(pmax)
plt.plot(w*1e-9,pflip1,'g-',label=r'$\omega_1=0.1GHz$')

pflip2 = []
w1 = 5e7
for i in w:   
    pmax = flip_prob([w0,i,w1],1.5e-7,750)
    pflip2.append(pmax)
#plot the datas
plt.plot(w*1e-9,pflip2,'r--',label=r'$\omega_1=0.05GHz$')
plt.xlim(0.9,1.1)
plt.ylim(0,1)
plt.xlabel(r'$\omega [GHz]$')
plt.ylabel(r'$P_{flip}$')
plt.title("Flipping Probability vs Frequency")
plt.legend()
plt.grid()
plt.show()
'''


def add_subplot_axes(ax,rect):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  
    #subax = fig.add_axes([x,y,width,height],facecolor=facecolor)  # matplotlib 2.0+
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

fig ,ax =plt.subplots()
N = 100_000
dt = 0.1e-9
t = np.linspace(0,N*dt,N)
y = integrate.odeint(up_down_prob,[1,0,0,0],t,args=([1e9,0.1e9,0.1e9],))
#print(y)
#ax.semilogy = ax.plot
pb1 = y[:,2]**2+y[:,3]**2
Yw = np.fft.fft(pb1)
Yw = np.abs(Yw)**2
wp = np.linspace(0,2*np.pi/dt,N)
ax.plot(wp*1e-9,Yw,label=r'$\omega=0.1GHz$')


y = integrate.odeint(up_down_prob,[1,0,0,0],t,args=([1e9,0.2e9,0.1e9],))
#print(y)
pb2 = y[:,2]**2+y[:,3]**2
Yw2 = np.fft.fft(pb2)
Yw2 = np.abs(Yw2)**2
ax.plot(wp*1e-9,Yw2,ls='--',label=r'$\omega=0.2GHz$')
ax.set_xlim(0,1.4)
ax.set_ylim(0,1.6e4)
ax.legend()


rect = [0.1,0.25,0.4,0.7]
subax = add_subplot_axes(ax,rect)
new_w = [0.6,0.7,0.8,0.9]
ls = ['-','--','-.',':']
for i,j in zip(new_w,ls):
    y = integrate.odeint(up_down_prob,[1,0,0,0],t,args=([1e9,i*1e9,0.1e9],))
    pb1 = y[:,2]**2+y[:,3]**2
    Yw = np.fft.fft(pb1)
    Yw = np.abs(Yw)**2
    subax.semilogy(wp*1e-9,Yw,ls=j,label=r'$\omega=%0.1fGHz$'%i)

subax.set_xlim(0,0.5)
subax.set_ylim(0.1,1e11)
subax.grid()
subax.legend()
ax.set_xlabel(r"$\omega^{\prime}$ [GHz]")
ax.set_ylabel("Power Spectrum (arb. units)")
plt.show()




