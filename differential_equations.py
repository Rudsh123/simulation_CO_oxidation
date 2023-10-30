import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def odes(x, t):
  Temp=539.5 #[K]
  R=8.31446261815324 #[J/mol*K]
  p_co=3 * 10**(-5) #[mbar]
  p_o=6.67 * 10**(-5) #[mbar]
  s_c=1
  k_c=3.135 * 10**(5) #[s^-1*mbar^-1]
  k_o=5.858 * 10**(5) #[s^-1*mbar^-1]
  u_s=1
  v_s=0.8
  k_r=3 * 10**(6) * np.exp(- ((10*4184)/(R*Temp)))  # [s^-1]
  k_d=2 * 10**(16) * np.exp(- ((38*4184)/(R*Temp))) # [s^-1]
  q=3

  so1=0.6
  so2=0.4
  kp=10**2 * np.exp(-7*4184/(Temp*R))
  r3=-(1/0.0135)
  r=[-0.026*r3, 0.3*r3, -1.05*r3, r3 ] #r0, r1, r2, r3
  poly_coefficent=0

  u = x[0]  # represents CO
  v = x[1]  # represents oxygen
  w = x[2]  # for the oscilating part

  if u<=0.2:
    dwdt=-w*kp
  elif u>=0.5:
    dwdt=(1-w)*kp
  else:
    for i in range(len(r)):
      poly_coefficent += r[i]*u**i
    dwdt = (poly_coefficent - w)*kp

  s_o = w*so1+(1-w)*so2
  dudt=p_co*k_c*s_c*(1-(u/u_s)**(q))-k_d*u-k_r*u*v
  dvdt=p_o*k_o*s_o*(1-u/u_s - v/v_s)**(2)-k_r*u*v

  return [dudt,dvdt,dwdt]


x0 = np.array([0,0,0])  # initial values (carbonyl/oxygen)
t = np.linspace(0, 20 , 100000)    # timeline

x=odeint(odes,x0,t)

u=x[:,0]
v=x[:,1]
w=x[:,2]

plt.plot(t, u, "b-")
plt.plot(t, v, "r-")
plt.plot (t , w,  "g--")
#plt.xscale("log")
plt.ylabel("blue= CO; red= O ; green =w")

plt.xlabel("Time [s]")
plt.show()
