# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 10:20:51 2023

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
# 1. E(w) with GDD0 
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


plt.subplot(212)
plt.plot(om[f1:f2], phi_om2[f1:f2], "b")
plt.xlabel("Angular frequency (PHz)")
plt.ylabel("Spectral phase (rad)")
plt.grid()
plt.tight_layout()

# =============================================================================
# I (w) and GD(w) of the output pulse 
# =============================================================================

I=np.abs(E_out)**2 
GD=np.gradient(phi_om2, dOM)

om1=om[I>max(I)/2][0]
om2=om[I>max(I)/2][-1]
# om_c0=om[I<=max(I)]
i1=np.argmin(om<=om1)
i2=np.argmin(om<=om2)
print("w1={}, w2={}".format(om1, om2))
print("i1={}, i2={}".format(i1, i2))

plt.figure(2)
plt.subplot(211)
plt.plot(om[f1:f2], I[f1:f2], "r")
plt.vlines(om1, 0, max(I)/2, color="red", linestyles="dotted")
plt.vlines(om2, 0, max(I)/2,  color="red", linestyles="dotted")
plt.hlines(max(I)/2, 1, om2, color="red", linestyles="dotted")
plt.ylabel("Spectral intensity (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[f1:f2], GD[f1:f2], "b")
plt.vlines(om1, -0.5, 1, color="red", linestyles="dotted")
plt.vlines(om2, -0.5, 1,  color="red", linestyles="dotted")
plt.hlines(GD[i1], 1, om2,  color="red", linestyles="dotted")
plt.hlines(GD[i2], 1, om2,  color="red", linestyles="dotted")
plt.xlabel("Angular frequency (PHz)")
plt.ylabel("Rel group delay (fs)")
plt.grid()
plt.tight_layout()

delta_GD=np.abs(GD[i1]-GD[i2])
print("delta GD={}".format(delta_GD))

# =============================================================================
# E(t) and phi(t)
# =============================================================================

E_0OM=np.fft.fftshift(E_out)
C_OT=(OM/(2*np.pi))*np.fft.ifft(E_0OM)
t_OT=np.fft.fftfreq(len(E_0OM), dOM/(2*np.pi))
C=np.fft.fftshift(C_OT)
t=np.fft.fftshift(t_OT)

#unwrap
phi_t=np.unwrap(np.angle(C))

i0=np.argmin(t<=-50)
i1=np.argmin(t<=50)

f0=np.argmin(np.abs(t))

n=int(phi_t[f0]/(2*np.pi))

phi_t2=phi_t-n*2*np.pi


plt.figure(3)
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
# I(t) and w_ins(t)
# =============================================================================

I_t=np.abs(C)**2 
w_ins=np.gradient(phi_t2, dt)

t1=t[I_t>=max(I_t)/2][0]
t2=t[I_t>=max(I_t)/2][-1]

plt.figure(4)
plt.subplot(211)
plt.plot(t[i0:i1], I_t[i0:i1], "r")
plt.vlines(t1, 0, max(I_t)/2, color="red", linestyles="dotted")
plt.vlines(t2, 0, max(I_t)/2,  color="red", linestyles="dotted")
plt.hlines(max(I_t)/2, t[i0], t2, color="red", linestyles="dotted")
plt.ylabel("Electric field (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(t[i0:i1], w_ins[i0:i1], "b")
plt.vlines(t1, w_ins[i0], w_ins[i1], color="red", linestyles="dotted")
plt.vlines(t2, w_ins[i0], w_ins[i1],  color="red", linestyles="dotted")
plt.ylabel("Instantaneous frequency (rad)")
plt.xlabel("Time (fs)")
plt.grid()

delta_t=np.abs(t1-t2)
anal_t=tau0*np.sqrt(1+(4*np.log(2)*GDD0/tau0**2)**2)

print("pulse duration from graph \n t={} \n analytical t={}".format(delta_t, anal_t))


# =============================================================================
# 1. E(w) with GDD0 =8000
# =============================================================================

i_min=np.argmax(om>=1)
i_max=np.argmin(om<=4)

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


plt.figure(1)
plt.subplot(211)
plt.plot(om[i_min:i_max], np.abs(E_out)[i_min:i_max], "r")
plt.ylabel("Spectral amplitude (a.u)")
plt.grid()

# =============================================================================
# I (w) and GD(w) of the output pulse 
# =============================================================================

I=np.abs(E_out)**2 
GD=np.gradient(phi_om2, dOM)

om1=om[I>max(I)/2][0]
om2=om[I>max(I)/2][-1]
# om_c0=om[I<=max(I)]
i1=np.argmin(om<=om1)
i2=np.argmin(om<=om2)
c0=n
print("w1={}, w2={}".format(om1, om2))
print("i1={}, i2={}".format(i1, i2))

plt.figure(2)
plt.subplot(211)
plt.plot(om[i_min:i_max], I[i_min:i_max], "r")
plt.vlines(om1, 0, max(I)/2, color="red", linestyles="dotted")
plt.vlines(om2, 0, max(I)/2,  color="red", linestyles="dotted")
plt.hlines(max(I)/2, 1, om2, color="red", linestyles="dotted")
plt.ylabel("Spectral intensity (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(om[i_min:i_max], GD[i_min:i_max], "b")
plt.vlines(om1, -0.5, 1, color="red", linestyles="dotted")
plt.vlines(om2, -0.5, 1,  color="red", linestyles="dotted")
plt.hlines(GD[i1], 1, om2,  color="red", linestyles="dotted")
plt.hlines(GD[i2], 1, om2,  color="red", linestyles="dotted")
plt.xlabel("Angular frequency (PHz)")
plt.ylabel("Rel group delay (fs)")
plt.grid()
plt.tight_layout()

delta_GD=np.abs(GD[i1]-GD[i2])
print("delta GD={}".format(delta_GD))

# =============================================================================
# E(t) and phi(t)
# =============================================================================

E_0OM=np.fft.fftshift(E_out)
C_OT=(OM/(2*np.pi))*np.fft.ifft(E_0OM)
t_OT=np.fft.fftfreq(len(E_0OM), dOM/(2*np.pi))
C=np.fft.fftshift(C_OT)
t=np.fft.fftshift(t_OT)

#unwrap
phi_t=np.unwrap(np.angle(C))

i0=np.argmin(t<=-50)
i1=np.argmin(t<=50)

f0=np.argmin(np.abs(t))

n=int(phi_t[f0]/(2*np.pi))

phi_t2=phi_t-n*2*np.pi


plt.figure(3)
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
# I(t) and w_ins(t)
# =============================================================================

I_t=np.abs(C)**2 
w_ins=np.gradient(phi_t2, dt)

t1=t[I_t>=max(I_t)/2][0]
t2=t[I_t>=max(I_t)/2][-1]

plt.figure(4)
plt.subplot(211)
plt.plot(t[i0:i1], I_t[i0:i1], "r")
plt.vlines(t1, 0, max(I_t)/2, color="red", linestyles="dotted")
plt.vlines(t2, 0, max(I_t)/2,  color="red", linestyles="dotted")
plt.hlines(max(I_t)/2, t[i0], t2, color="red", linestyles="dotted")
plt.ylabel("Electric field (a.u)")
plt.grid()

plt.subplot(212)
plt.plot(t[i0:i1], w_ins[i0:i1], "b")
plt.vlines(t1, w_ins[i0], w_ins[i1], color="red", linestyles="dotted")
plt.vlines(t2, w_ins[i0], w_ins[i1],  color="red", linestyles="dotted")
plt.ylabel("Instantaneous frequency (rad)")
plt.xlabel("Time (fs)")
plt.grid()

delta_t=np.abs(t1-t2)
anal_t=tau0*np.sqrt(1+(4*np.log(2)*GDD0/tau0**2)**2)

print("pulse duration from graph \n t={} \n analytical t={}".format(delta_t, anal_t))

