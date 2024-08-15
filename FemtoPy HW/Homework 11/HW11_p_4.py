# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 22:19:21 2023

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
# 9. GDD0=8000 
# =============================================================================

GDD0=8000

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

#Inverse fourier transform: 

E_OT=np.fft.fftshift(E_out)
E_t=(N/(2*np.pi))*np.fft.ifft(E_OT)
t_OT=np.fft.fftfreq(len(E_OT), dOM/(2*np.pi))
E_t=np.fft.fftshift(E_t)
t=np.fft.fftshift(t_OT)


#unwrap the angle: 

it0=np.argmin(abs(t))
    
tem_phase=np.unwrap(np.angle(-E_t))
n=int(tem_phase[it0]/(2*np.pi))
tem_phase2=tem_phase-n*np.pi

I_t=np.abs(E_t)**2 
w_ins=np.gradient(tem_phase, dt)

plt.figure(9)
plt.subplot(211)
plt.plot(t, I_t, "r")
plt.ylabel("Electric field (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(t, w_ins, "b")
plt.ylabel("Instantanous angular frequency")
plt.xlabel("Time, (fs)")
plt.grid()
plt.tight_layout()
plt.show()

# =============================================================================
# 10. Increasing the frequency of the array elements 
# =============================================================================


N=2**18
dOM=OM/N
om=np.arange(-25, 25, dOM)

E_in=(np.sqrt(np.pi)*tau0/(np.sqrt(8*np.log(2))))*np.exp(-(om-om0)**2*tau0**2/(8*np.log(2)))

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

#Inverse fourier transform: 

E_OT=np.fft.fftshift(E_out)
E_t=(N/(2*np.pi))*np.fft.ifft(E_OT)
t_OT=np.fft.fftfreq(len(E_OT), dOM/(2*np.pi))
E_t=np.fft.fftshift(E_t)
t=np.fft.fftshift(t_OT)


#unwrap the angle: 

it0=np.argmin(abs(t))
    
tem_phase=np.unwrap(np.angle(-E_t))
n=int(tem_phase[it0]/(2*np.pi))
tem_phase2=tem_phase-n*np.pi

I_t=np.abs(E_t)**2 
w_ins=np.gradient(tem_phase, dt)

plt.figure(10)
plt.subplot(211)
plt.plot(t, I_t, "r")
plt.ylabel("Temporal intensity (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(t, w_ins, "b")
plt.ylabel("Instantanous angular frequency")
plt.xlabel("Time, (fs)")
plt.grid()
plt.tight_layout()
plt.show()
