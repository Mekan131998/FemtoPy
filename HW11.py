# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 13:56:39 2023

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

E_in=(np.sqrt(np.pi)*tau0/(np.sqrt(8*np.log(2))))*np.exp(-(om-om0)**2*tau0**2/(8*np.log(2)))

phase=np.angle(E_in)

# =============================================================================
# 1. E(w) with GDD0=+50
# =============================================================================
GDD0=50

phi_w=0.5*GDD0*(om-om0)**2

h=np.exp(-1j*phi_w)

E_out=h*E_in
phase_out=-np.angle(E_out)

phi_om=np.unwrap(phase_out)

f1=np.argmin(np.abs(om-1))
f2=np.argmin(np.abs(om-4))

fr0=np.argmin(np.abs(om-om0))
n=int(phi_om[fr0]/(2*np.pi))

phi_om2=phi_om-n*2*np.pi


plt.figure(1)
plt.subplot(211)
plt.plot(om[f1:f2], np.abs(E_out)[f1:f2], "r")
plt.ylabel("Spectral amplitude (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[f1:f2], phi_om2[f1:f2], "b")
plt.xlabel("Angular frequency (PHz)")
plt.ylabel("Spectral phase (rad)")
plt.grid()
plt.tight_layout()

#################################################################################
#Function for finding FWHM values of intensity 
# arg_FWHM can be either time or angular frequency 
# it returns FWHM values of time or angular frequency and corresponding idexies  

def FWHM(intensity, arg_FWHM):
    t1, t2=arg_FWHM[intensity>max(intensity)/2][0], arg_FWHM[intensity>max(intensity)/2][-1]
    i1, i2=np.argmin(arg_FWHM<t1), np.argmin(arg_FWHM<t2)
    return t1, t2, i1, i2

##############################################################################


# =============================================================================
# 2. I (w) and GD(w) of the output pulse 
# =============================================================================

I=np.abs(E_out)**2 
GD=np.gradient(phi_om2, dOM)

#finding angular frequency values: 

w1, w2, i1, i2=FWHM(I, om)

plt.figure(2)
plt.subplot(211)
plt.plot(om[f1:f2], I[f1:f2], "r")
plt.vlines(w1, 0, max(I)/2, linestyles="dashed", color='red')
plt.vlines(w2, 0, max(I)/2, linestyles="dashed", color='red')
plt.ylabel("Spectral intensity (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[f1:f2], GD[f1:f2], "b")
plt.vlines(w1, min(GD[f1:f2]), max(GD[f1:f2]), linestyles="dashed", color='red')
plt.vlines(w2, min(GD[f1:f2]), max(GD[f1:f2]), linestyles="dashed", color='red')
plt.hlines(GD[i1], w1, w2, linestyles="dashed", color="blue")
plt.hlines(GD[i2], w1, w2, linestyles="dashed", color="blue")
plt.xlabel("Angular frequency, PHz")
plt.ylabel("Rel. group delay (fs)")
plt.grid()
plt.tight_layout()
plt.show()

d_GD=np.abs(GD[i1]-GD[i2])
print("intexes \t i1={}, \t i2={}".format(i1, i2))
print("\t w1={0:0.2f}, \t w2={1:0.2f}".format(w1, w2))
print("Delta GD={0:0.2f}".format(d_GD))

# =============================================================================
# 3. E(t) AND ϕ(t) OF THE OUTPUT PULSE
# =============================================================================

#Inverse fourier transform: 

E_OT=np.fft.fftshift(E_out)
E_t=(N/(2*np.pi))*np.fft.ifft(E_OT)
t_OT=np.fft.fftfreq(len(E_OT), dOM/(2*np.pi))
E_t=np.fft.fftshift(E_t)
t=np.fft.fftshift(t_OT)

it1=np.argmin(t<=-50)
it2=np.argmin(t<=50)

#unwrap the angle: 

it0=np.argmin(abs(t))
    
tem_phase=np.unwrap(np.angle(-E_t))
n=int(tem_phase[it0]/(2*np.pi))
tem_phase2=tem_phase-n*np.pi


plt.figure(3)
plt.subplot(211)
plt.plot(t[it1:it2], E_t[it1:it2], "r")
plt.plot(t[it1:it2], np.abs(E_t[it1:it2]), "r")
plt.ylabel("Electric field (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(t[it1:it2], tem_phase2[it1:it2], "b")
plt.ylabel("Temporal phase (a.u)")
plt.xlabel("Time, (fs)")
plt.grid()
plt.tight_layout()
plt.show()

# =============================================================================
#4. I(t) AND ωins (t) OF THE OUTPUT PULSE
# =============================================================================

I_t=np.abs(E_t)**2 
w_ins=np.gradient(tem_phase, dt)

#finding FWHM times: 
t1, t2, idx1, idx2=FWHM(I_t, t)

#analytical formula for the pulse duration: 
    
t_a=tau0*np.sqrt(1+((4*np.log(2)*GDD0)/tau0**2)**2)

plt.figure(4)
plt.subplot(211)
plt.plot(t[it1:it2], I_t[it1:it2], "r")
plt.ylabel("Electric field (a.u)")
plt.vlines(t1, 0, max(I_t)/2, linestyles="dashed", color='red')
plt.vlines(t2, 0, max(I_t)/2, linestyles="dashed", color='red')
plt.hlines(max(I_t)/2, t[it1], t2, linestyles="dashed", color='red')
plt.grid()

plt.subplot(212)
plt.plot(t[it1:it2], w_ins[it1:it2], "b")
plt.vlines(t1, min(w_ins[it1:it2]), max( w_ins[it1:it2]), linestyles="dashed", color='red')
plt.vlines(t2, min(w_ins[it1:it2]), max( w_ins[it1:it2]), linestyles="dashed", color='red')
plt.ylabel("Instantanous angular frequency")
plt.xlabel("Time, (fs)")
plt.grid()
plt.tight_layout()
plt.show()

d_t=np.abs(t1-t2)
print("intexes \t i1={}, \t i2={}".format(idx1, idx2))
print("\t t1={0:0.2f}, \t t2={1:0.2f}".format(t1, t2))
print("Delta \t t={0:0.2f}".format(d_t))
print("Analytical \t t={0:0.2f}".format(t_a))


