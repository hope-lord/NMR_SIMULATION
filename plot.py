from numpy import abs
from nmr import  prob_down, power_spectrum, diff_eq_of_spin
import matplotlib.pylab as plt

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

######### UPDOWN PROBABILITY
w0 = 1e9
w = 1e9
w1 =0.1e9
t , y = prob_down(0.1e-9,2000,w,w1,w0)
plt.plot(plt.array(t)*1e9,y,'r-')
plt.xlabel('Time [ns]')
plt.ylabel(r'$P_b(t)$')
plt.xlim(0,200)
plt.ylim(0,1)
plt.show()
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

fig ,ax =plt.subplots()
N = 2**17  #2^17
dt = 0.1e-9
wp,Yw =  power_spectrum(dt,N,0.1e9,0.1e9,1e9)

#ax.semilogy = ax.plot

ax.plot(plt.array(wp)*1e-9,Yw,label=r'$\omega=0.1GHz$')

wp,Yw2 =  power_spectrum(dt,N,0.2e9,0.1e9,1e9)

ax.plot(plt.array(wp)*1e-9,Yw2,ls='--',label=r'$\omega=0.2GHz$')
ax.set_xlim(0,1.4)
ax.set_ylim(0,5e4)
ax.legend()


rect = [0.1,0.25,0.4,0.7]
subax = add_subplot_axes(ax,rect)
new_w = [0.6,0.7,0.8,0.9]
ls = ['-','--','-.',':']
for i,j in zip(new_w,ls):
    wp,Yw = power_spectrum(dt,N,i*1e9,0.1e9,1e9)
    #Yw = abs(plt.fft(Yw))**2
    subax.semilogy(plt.array(wp)*1e-9,Yw,ls=j,label=r'$\omega=%0.1fGHz$'%i)

subax.set_xlim(0,0.5)
subax.set_ylim(0.1,1e11)
subax.grid()
subax.legend()
ax.set_xlabel(r"$\omega^{\prime}$ [GHz]")
ax.set_ylabel("Power Spectrum (arb. units)")
plt.show()




