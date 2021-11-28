from math import pi,cos, sin

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
        return part1+part2
