# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
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

# =============================================================================
# 6. Calculation the effect of GD_0
# E_out_GD (w) if phi(w)=GD_0*(w-w0)
# =============================================================================
omega0=2.3546
GD_0=20 
phi=GD_0*(om-omega0)

H=np.exp(-1j*phi)

E_out_GD=H*E_in 

phi_om=np.angle(E_out_GD)

#finding indexes close to the 1Hz and 4Hz 
phi_om=np.unwrap(-phi_om)

f1=np.argmin(np.abs(om-1))
f2=np.argmin(np.abs(om-4))

fr0=np.argmin(np.abs(om-omega0))
n=int(phi_om[fr0]/(2*np.pi))

phi_om2=phi_om-n*2*np.pi

plt.figure(6)
plt.subplot(211)
plt.plot(om[f1:f2], np.abs(E_out_GD)[f1:f2], "r")
plt.ylabel("Amplitude (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[f1:f2], phi_om2[f1:f2], "b")
plt.xlabel("angular frequency (PHz)")
plt.ylabel("Spectral phase (rad) ")
plt.grid()
plt.tight_layout()

# =============================================================================
# 7. I (w) and GD(w) of the output pulse
# =============================================================================

I_w=np.abs(E_out_GD)**2

GD_w=np.gradient(phi_om2, dOM)

plt.figure(7)
plt.subplot(211)
plt.plot(om[f1:f2], I_w[f1:f2], "r")
plt.ylabel("Spectral intensity (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[f1:f2], GD_w[f1:f2], "b")
plt.xlabel("Angular frequency (PHz) ")
plt.ylabel("Rel. group delay (fs)")
plt.grid()
plt.tight_layout()

# =============================================================================
# 8. E (t) and phi (t) of the output pulse 
# =============================================================================

E_0OM=np.fft.fftshift(E_out_GD)
C_OT=(OM/(2*np.pi))*np.fft.ifft(E_0OM)
t_OT=np.fft.fftfreq(len(E_0OM), dOM/(2*np.pi))
C=np.fft.fftshift(C_OT)
t=np.fft.fftshift(t_OT)

#unwrap
phi_t=np.unwrap(np.angle(C))

i0=np.argmin(t<=10)
i1=np.argmin(t<=30)

f0=np.argmin(np.abs(t))

n=int(phi_t[f0]/(2*np.pi))

phi_t2=phi_t-n*2*np.pi


plt.figure(8)
plt.subplot(211)
plt.plot(t[i0:i1], C[i0:i1], "r")
plt.plot(t[i0:i1], np.abs(C)[i0:i1], "r--")
plt.ylabel("Electric field (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(t[i0:i1], phi_t2[i0:i1], "b")
plt.ylabel("Phase (rad)")
plt.xlabel("Time (fs)")
plt.grid()

# =============================================================================
#9. Intensity I(t) and instantaneous angular frequency w_in(t) as a function of time
# =============================================================================

I_t=np.abs(C)**2 
w_ins=np.gradient(phi_t2, dt)


plt.figure(9)
plt.subplot(211)
plt.plot(t[i0:i1], I_t[i0:i1], "r")
plt.ylabel("Intensity (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(t[i0:i1], w_ins[i0:i1], "b")
plt.ylabel("Ins. ang. freq. (PHz)")
plt.xlabel("Time (fs)")
plt.grid()
plt.tight_layout()

