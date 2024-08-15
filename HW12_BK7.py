# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 22:13:46 2023

@author: Owner
"""

import numpy as np 
import matplotlib.pyplot as plt 
plt.close('all')

c=299.792  #nm/fs
lmbda=800 
tau0=6 
om0=2*np.pi*c/lmbda

OM=50
N=2**16
dOM=OM/N
om=np.arange(-OM/2, OM/2, dOM)
dt=2*np.pi/OM 

# =============================================================================
# SF10 gType=1,  range: 0.38, 2.5
# BK7 gType=2, range:0.3, 2.5
# FS gType=3, range:0.21, 6.7
# =============================================================================

def refindex(an_freq, gType):
    if an_freq==0:
        x=0.8           #um
    else:
        x=(2*np.pi*c/an_freq)*1e-3
    if gType==1 and x>0.38 and x<2.5:
        return (1+1.62153902/(1-0.0122241457/x**2)+0.256287842/(1-0.0595736775/x**2)+1.64447552/(1-147.468793/x**2))**.5
    elif gType==2 and x>0.3 and x<2.5:
        return (1+1.03961212/(1-0.00600069867/x**2)+0.231792344/(1-0.0200179144/x**2)+1.01046945/(1-103.560653/x**2))**.5    
    elif gType==3 and x>0.21 and x<6.7:
        return (1+0.6961663/(1-(0.0684043/x)**2)+0.4079426/(1-(0.1162414/x)**2)+0.8974794/(1-(9.896161/x)**2))**.5
    else: return 1

#Test function: 
def testing(lmda, gType):   #here lmbda in mkm (1e06 m)
    lmbda1=lmda*1e3      #here lmda in nm
    an_freq1=2*np.pi*c/lmbda1      # an frequency in PHz
    n1=refindex(an_freq1, gType)
    print("calculates \t", n1)

print("Should return \t 1.7113")
testing(0.8, 1)
print("Should return \t 1.5131")
testing(0.7, 2)
print("Should return \t  1.4504")
testing(1, 3)
print("Should return 1")
testing(0.1, 1)
testing(3, 2)
testing(7, 3)

#Optical transfer function for BK7 gType=2: 
L=1*1e6        #in nm 
ref_FS=[]

for i in om: 
    ref_FS.append(refindex(i, 2))
    
H_w=np.exp(-1j*(om/c)*ref_FS*L)

E_in=np.sqrt(np.pi)*tau0/np.sqrt(8*np.log(2))*np.exp(-(om-om0)**2*tau0**2/(8*np.log(2)))

E_out=H_w*E_in

phase_out=-np.angle(E_out)
phi_om=np.unwrap(phase_out)
fr0=np.argmin(np.abs(om-om0))
n=int(phi_om[fr0]/(2*np.pi))
phi_om2=phi_om-n*2*np.pi

i1=np.argmin(np.abs(om-1))
i2=np.argmin(np.abs(om-4))


plt.figure()
plt.subplot(211)
plt.plot(om[i1:i2], abs(E_out)[i1:i2], "r")
plt.ylabel("Spectral amplitude (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[i1:i2], -phi_om2[i1:i2], "b")
plt.ylabel("Spectral phase (rad)")
plt.xlabel("Angular frequency, PHz")
plt.grid()
plt.tight_layout()
plt.show()

#intensity:
I=np.abs(E_out)**2 

#group delay: 
GD=np.gradient(phi_om2, dOM)

def FWHM(intensity, arg_FWHM):
    t1, t2=arg_FWHM[intensity>max(intensity)/2][0], arg_FWHM[intensity>max(intensity)/2][-1]
    i1, i2=np.argmin(arg_FWHM<t1), np.argmin(arg_FWHM<t2)
    return t1, t2, i1, i2

w1, w2, h1, h2=FWHM(I, om)

print("indicies for FWHM: \n h1={}, \t h2={}".format(h1, h2))

delta_GD=np.abs(GD[h1]-GD[h2])
print('Delta GD:', delta_GD)

print("value of GD which is closest to w0: {}".format(GD[fr0]))

rel_GD=GD-GD[fr0]

plt.figure()
plt.subplot(211)
plt.plot(om[i1:i2], I[i1:i2], "r")
plt.vlines(om[h1], 0, max(I)/2, linestyles="dashed", color="red")
plt.vlines(om[h2], 0, max(I)/2, linestyles="dashed", color="red")
plt.ylabel("Spectral intensity (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[i1:i2], rel_GD[i1:i2])
plt.ylabel("Rel. group delay")
plt.xlabel("Angular frequency, (PHz)")
plt.grid()
plt.tight_layout()
plt.show()


GDD = np.gradient(GD, dOM) # Group delay dispersion [fs^2]
GDD0 = GDD[fr0]
print('Fused silica GDD(om0)= {0:0.2f} fs^2'.format(GDD0))
pulse_duration= tau0 * np.sqrt(1 + (4 * np.log(2) * GDD0/(tau0 ** 2))**2)
print('Pulse duration = {0:0.3f} fs'.format(pulse_duration))
# =============================================================================
# E(t) and phi(t) of the output pulse 
# =============================================================================

E_0_OM=np.fft.fftshift(E_out)
E_t0=(OM/(2*np.pi))*np.fft.ifft(E_0_OM)
t_0=np.fft.fftfreq(len(E_0_OM), dOM/(2*np.pi))
E_t=np.fft.fftshift(E_t0)
t=np.fft.fftshift(t_0)

tmax = np.argmin(t<(GD[fr0] + 2 * pulse_duration))
tmin = np.argmin(t<(GD[fr0] - 2 * pulse_duration))

#unwrap:
phi_om=np.unwrap(np.angle(E_t))
t00= np.argmin(t<GD[fr0])
n=int(phi_om[t00]/(2*np.pi))
phi_om2=phi_om-n*2*np.pi

plt.figure()
plt.subplot(211)
plt.plot(t[tmin:tmax], np.abs(E_t[tmin:tmax]), "r--")
plt.plot(t[tmin:tmax], E_t[tmin:tmax], "r")
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

plt.figure()
plt.subplot(211)
plt.plot(t[tmin:tmax], I_t[tmin:tmax], "r")
plt.vlines(half1, 0, max(I_t)/2, linestyles="dashed", color="red")
plt.vlines(half2, 0, max(I_t)/2, linestyles="dashed", color="red")
plt.hlines(max(I_t)/2, t[idx1], t[idx2], linestyles="dashed", color="red")
plt.ylabel("Temporal intensity")
plt.grid()


plt.subplot(212)
plt.plot(t[tmin:tmax], w_ins[tmin:tmax])
plt.vlines(half1, min(w_ins[tmin:tmax]), max(w_ins[tmin:tmax]), linestyles="dashed", color="red")
plt.vlines(half2, min(w_ins[tmin:tmax]), max(w_ins[tmin:tmax]), linestyles="dashed", color="red")
plt.ylabel("Instantaneous frequency")
plt.xlabel("Time, (fs)")
plt.grid()
plt.show()

delta_t=np.abs(half1-half2)

print("Duration from intensity: {0:0.3f}".format(delta_t))
