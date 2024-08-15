# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 11:08:41 2023

@author: Owner
"""


import numpy as np 
import matplotlib.pyplot as plt 

plt.close('all')
lmbda=800 
tau0=6 

OM=50
N=2**16
dOM=OM/N

om=np.arange(-25, 25, dOM)
dt=2*np.pi/OM 

c=299.73 

om0=2*np.pi*c/lmbda
T=2*np.pi/OM
t=np.arange(-T/2, T/2, dt)
E_in=np.exp(-(om-om0)**2*tau0**2/(8*np.log(2)))
L=1*1e10        #in nm 
ref_FS=[]

def refindex(an_freq):
    if an_freq==0:
        x=0.8           #um
    else:
        x=(2*np.pi*c/an_freq)*1e-3
    if x>0.23 and x<1.69:
        return 1+0.05792105/(238.0185-x**(-2))+0.00167917/(57.362-x**(-2))
    else: return 1

for i in om: 
    ref_FS.append(refindex(i))
   
H_w=np.exp(-1j*(om/c)*ref_FS*L)
E_out=H_w*E_in

phi_om0=-np.angle(E_out)
phi_om=np.unwrap(phi_om0)
fr0=np.argmin(np.abs(om-om0))
n=int(phi_om[fr0]/(2*np.pi))

phi_om2=phi_om-n*2*np.pi

GD=np.gradient(phi_om2, dOM)


def FWHM(intensity, arg_FWHM):
    t1, t2=arg_FWHM[intensity>max(intensity)/2][0], arg_FWHM[intensity>max(intensity)/2][-1]
    i1, i2=np.argmin(arg_FWHM<t1), np.argmin(arg_FWHM<t2)
    return t1, t2, i1, i2

print("w0 index {}".format(fr0))
print("GD {}".format(GD[fr0]))

phi_mod=phi_om2-GD[fr0]*(om-om0)


GD_r=np.gradient(phi_mod, dOM)

f1=np.argmin(np.abs(om-1.5))
f2=np.argmin(np.abs(om-3.25))

I=np.abs(E_in)**2

w1=om0-2*np.log(2)/tau0
w2=om0+2*np.log(2)/tau0

print("w1={}, w2={}".format(w1, w2))

i1=np.argmin(om<w1)
i2=np.argmin(om<w2)


plt.figure(1)
plt.subplot(211)
plt.plot(om[f1:f2], I[f1:f2], "r")
plt.vlines(w1, 0, max(I)/2, linestyles="dashed", color='red')
plt.vlines(w2, 0, max(I)/2, linestyles="dashed", color='red')
plt.ylabel("Spectral intensity (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[f1:f2], GD_r[f1:f2], "b")
plt.xlabel("Angular frequency, PHz")
plt.ylabel("Rel. group delay (fs)")
plt.grid()
plt.tight_layout()
plt.show()

print("delta GD={}".format(GD[i2]-GD[i1]))



E_0_OM=np.fft.fftshift(E_out)
E_t0=(OM/(2*np.pi))*np.fft.ifft(E_0_OM)
t_0=np.fft.fftfreq(len(E_0_OM), dOM/(2*np.pi))
t=np.fft.fftshift(t_0)
E_t=np.fft.fftshift(E_t0)


tmin=np.argmin(t<-200)
tmax=np.argmin(t<200)

#unwrap:
phi_om=np.unwrap(-np.angle(E_t))
n=int(phi_om[fr0]/(2*np.pi))
phi_om2=phi_om-n*2*np.pi

I_t=np.abs(E_t)**2

plt.figure(2)
plt.subplot(211)
plt.plot(t[tmin:tmax], I_t[tmin:tmax], "r--")
plt.ylabel("Electric field (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(t[tmin:tmax], phi_om2[tmin:tmax])
plt.ylabel("Temporal phase (rad)")
plt.xlabel("Time (fs)")
plt.grid()
plt.tight_layout()
plt.show()

#intensity: 
I_t=np.abs(E_t)**2 
w_ins=np.gradient(phi_om2, dt)

half1, half2, idx1, idx2=FWHM(I_t, t)

plt.figure(3)
plt.subplot(211)
plt.plot(t[tmin:tmax], I_t[tmin:tmax], "r")
plt.ylabel("Temporal intensity")
plt.grid()

plt.subplot(212)
plt.plot(t[tmin:tmax], w_ins[tmin:tmax])
plt.ylabel("Instantaneous frequency")
plt.xlabel("Time, (fs)")
plt.grid()
plt.show()

delta_t=np.abs(half1-half2)

print("Duration from intensity: {0:0.3f}".format(delta_t))


