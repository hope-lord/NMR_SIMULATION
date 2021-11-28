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