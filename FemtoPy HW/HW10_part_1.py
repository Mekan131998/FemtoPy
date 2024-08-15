# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 06:38:41 2023

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

c=299.73 

om0=2*np.pi*c/lmbda

E_in=(np.sqrt(np.pi)*tau0/(np.sqrt(8*np.log(2))))*np.exp(-(om-om0)**2*tau0**2/(8*np.log(2)))

phase=np.angle(E_in)

# =============================================================================
# 1. E(w) of the input pulse 
# =============================================================================

i_min=np.argmax(om>=1)
i_max=np.argmin(om<=4)

plt.figure(1)
plt.subplot(211)
plt.plot(om[i_min:i_max], np.abs(E_in)[i_min:i_max], "r")
plt.ylabel("Spectral amplitude (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[i_min:i_max], phase[i_min:i_max], "b")
plt.xlabel("Angular frequency (pHz)")
plt.ylabel("Spectral phase (rad)")
plt.grid()
plt.tight_layout()

# =============================================================================
# 2. E(w) of the output pulse 
# =============================================================================

phi0=2
h=np.exp(-1j*phi0)

E_out=h*E_in
phase_out=-np.angle(E_out)

plt.figure(2)
plt.subplot(211)
plt.plot(om[i_min:i_max], np.abs(E_out)[i_min:i_max], "r")
plt.ylabel("Spectral amplitude (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[i_min:i_max], phase_out[i_min:i_max], "b")
plt.xlabel("Angular frequency (PHz)")
plt.ylabel("Spectral phase (rad)")
plt.grid()
plt.tight_layout()

# =============================================================================
# 3. I(w) and GD(w) of the output pulse
# =============================================================================
GD=np.gradient(phase_out, dOM)
IW=np.abs(E_out)**2
plt.figure(3)
plt.subplot(211)
plt.plot(om[i_min:i_max], IW[i_min:i_max], "r")
plt.ylabel("Spectral intenstiy (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[i_min:i_max], GD[i_min:i_max], "b")
plt.xlabel("Angular frequency (pHz)")
plt.ylabel("Rel. group delay (fs)")
plt.grid()
plt.tight_layout()

# =============================================================================
# 4. E(t) and phi(t) of the output pulse
# =============================================================================

# 4.1. Inverse Fourier transform

E_0OM=np.fft.fftshift(E_out)
C_0OM=(OM/(2*np.pi))*np.fft.ifft(E_0OM)
t_0OM=np.fft.fftfreq(len(E_0OM), dOM/(2*np.pi))
C=np.fft.fftshift(C_0OM)
t=np.fft.fftshift(t_0OM)

#4.2 Unwrap

temp_phase=np.unwrap(np.angle(C))

i0=np.argmin(t<=-15)
i1=np.argmin(t<=15)

f0=np.argmin(np.abs(t))

n=int(temp_phase[f0]/(2*np.pi))

temp_phase2=temp_phase-n*2*np.pi


plt.figure(4)
plt.subplot(211)
plt.plot(t[i0:i1], C[i0:i1], "r")
plt.plot(t[i0:i1], np.abs(C)[i0:i1], "r--")
plt.ylabel("Electric field (a.u)")
plt.xlim(-15,15)
plt.grid()

plt.subplot(212)
plt.plot(t[i0:i1], temp_phase2[i0:i1], "b")
plt.xlabel("Time (fs)")
plt.ylabel("Phase (rad)")
plt.xlim(-15, 15)
plt.grid()
plt.tight_layout()

#time step 
dt=2*np.pi/OM 

# =============================================================================
# 5. Intensity I(t) and the instantaneous w_ins(t) of the output pulse 
# =============================================================================

I=np.abs(C)**2
om_ins=np.gradient(temp_phase2, dt)

plt.figure(5)
plt.subplot(211)
plt.plot(t[i0:i1], I[i0:i1], "r")
plt.ylabel("Intensity (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(t[i0:i1], om_ins[i0:i1], "b")
plt.xlabel("Time (fs) ")
plt.ylabel("Instantaneous ang. frequency (PHz)")
plt.grid()
plt.tight_layout()
