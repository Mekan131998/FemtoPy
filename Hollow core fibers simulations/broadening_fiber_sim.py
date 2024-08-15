# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:16:50 2024

@author: Mekan Merdanov 

Simulation of the broadening through fiber 

using article: R. H. Stolen, C. Lin, “Self-phase-modulation in silica optical fibers,” Phys. Rev. A,
vol. 17, pp. 1448–1453, 1978

book: Ursula Keller Ultrafast lasers
"""

import numpy as np 
import matplotlib.pyplot as plt 

#%% Parameters and formulas: 
tau0=180e-12      # s, pulse duration 

T=10000e-12         # s, time period 
N=2**12      # number of points

dt=T/N      # number of steps 

t=np.arange(-T/2, T/2, dt)   # time array 

P_avg=90*1e-3         # W, average power for phi_max=4.5 pi

duty_cycle=90*1e-3/3
P_peak=P_avg/duty_cycle        # peak power

P_t=P_peak*np.exp(-4*np.log(2)*(t/tau0)**2)        # Gaussian pulse 

c=3*1e8;          # m/s, speed of ligh 
A_eff=9.61*1e-12   # m^2 effective area 
L_eff=81.2        # m, effective length
n=1.2             # refractive index

E2=8*np.pi*P_t/(n*c*A_eff);

n2_esu=1.16*1e-13     # esu, coefficient of self focusing
n2=4.19*1e-7*n2_esu/n # m^2/W

dn=(5/6)*0.5*n2*E2  # intensity dependent refractive index 

lmbda=514.5*1e-9    # m, central wavelength
om0=2*np.pi*c/lmbda # s^-1, circular frequency 

T0=2*np.pi/om0 

dphi_t=2*np.pi*L_eff*dn/lmbda;


U_t=np.sqrt(P_t)*np.exp(1j*(dphi_t)*t);  # Function that is needed to be taken FFT 

#%% Phase plot 

plt.figure(1)
plt.plot(t, dphi_t)
plt.show()

om_t=om0+np.gradient(dphi_t, t)

plt.figure(2)
plt.plot(t, om_t)
plt.show()


#%% Performing Fourier Transform

Et_0T=np.fft.fftshift(U_t)
C_0OM=np.fft.fft(Et_0T)/N

freq_0OM=np.fft.fftfreq(len(Et_0T), dt)

C=np.fft.fftshift(C_0OM)
freq=np.fft.fftshift(freq_0OM);

om=2*np.pi*freq


#%% Plotting the spectrum: 
    
plt.figure(3)
plt.plot(om, abs(C))
plt.show()














