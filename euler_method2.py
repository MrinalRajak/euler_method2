
# euler method2
# Solution of second order differential equation in RK4 method and
# comparison of same using scipy.integrate,odeint
# Test problem Damped Harmonic oscillator m(d/dt(dx/dt))+ lam(dx/dt)+kx=0
# we can write (dx/dt)=v and (dv/dt)=(-lam*v-k*x)/m
# initial condition : at t=0,x=0,v=1

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

x=0.0
v=1.0
I=np.array([x,v])
ti=0.0
tf=20.0
h=0.001
t=ti
k=7
lam=0.5
m=1.0
ts=[];xs=[];vs=[]
def Rk4(I,t):
    k1=dI(I,t)
    k2=dI(I+h/2.0*k1,t+h/2.0)
    k3=dI(I+h/2.0*k2,t+h/2.0)
    k4=dI(I+h*k3,t+h)
    I+=(k1+2*k2+2*k3+k4)/6.0*h
    t+=h
    return I,t
def dI(I,t):
    I=x,v
    dxdt=v
    dvdt=(-lam*v-k*x)/m
    return np.array([dxdt,dvdt])

while(t<tf):
    I,t=Rk4(I,t)
    x,v=I
    ts.append(t)
    xs.append(x)
    vs.append(v)

#solving the problem using scipy.integrate odeint
u0=[0.0,1] # initial value of and dxdt
tode=np.arange(ti,tf,h)
def damped_oscillator(u,tode):
    x1=u[0]
    v1=u[1]
    dxdt=v1
    dvdt=(-lam*v1-k*x1)/m
    return np.array([dxdt,dvdt])
sol=odeint(damped_oscillator,u0,tode)
x1=sol[:,0]
v1=sol[:,1]
plt.plot(ts,xs)
plt.plot(ts,vs)
plt.plot(tode,x1,label='velocity',ls='--',lw=2,c='red')
plt.plot(tode,v1,label='acceleration',ls='--',lw=2,c='blue')
plt.legend()
plt.grid(True)
plt.show()


















