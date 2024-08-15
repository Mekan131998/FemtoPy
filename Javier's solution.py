# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 23:13:53 2023

@author: Owner
"""

import numpy as np
import matplotlib.pyplot as plt

def refrindex(w, glass):
    if w == 0:
        x = 1
    else:
        x = 2 * np.pi * c / (w * 1000) # Wavelength [um]
    n = 1
    if glass == 'FS' and x > 0.24 and x < 6.7:
        n=(1+0.6961663/(1-(0.0684043/x)**2)+0.4079426/(1-(0.1162414/x)**2)+0.8974794/(1-(9.896161/x)**2))**.5
    elif glass == 'BK7' and x > 0.3 and x < 2.5:
        n=(1+1.03961212/(1-0.00600069867/x**2)+0.231792344/(1-0.0200179144/x**2)+1.01046945/(1-103.560653/x**2))**.5
    elif glass == 'SF10' and x > 0.38 and x < 2.5:
        n=(1+1.62153902/(1-0.0122241457/x**2)+0.256287842/(1-0.0595736775/x**2)+1.64447552/(1-147.468793/x**2))**.5
    return n

wl0 = 800
c = 299.792458
omega = 50
N = 2**17
dw = omega / N
w = np.arange(-omega/2, omega/2, dw)
w0 = 2 * np.pi * c / wl0
tau_0 = 6
T = 2 * np.pi / dw

n_800_fused_silica = refrindex(w0, 'fused silica')
n_800_BK7 = refrindex(w0, 'BK7')
n_800_SF10 = refrindex(w0, 'SF10')
print('n fused silica at 800 nm = ', n_800_fused_silica)
print('n BK7 at 800 nm = ', n_800_BK7)
print('n SF10 at 800 nm = ', n_800_SF10)


Ew = np.sqrt(np.pi) * tau_0 / (np.sqrt(8 * np.log(2))) * np.exp(-(w - w0)**2 * tau_0**2 / (8 * np.log(2)))
L = 1e6 # [nm]
n_w_FS = [refrindex(i, 'FS') for i in w]
E_out = np.exp(-1j * (w/c) * n_w_FS * L) * Ew

phase = -np.angle(E_out)
spectral_phase = np.unwrap(phase)
indexw0 = np.argmin(np.abs(w-w0))
n = int(spectral_phase[indexw0]/(2 * np.pi))
phase_shifted = spectral_phase  - n * 2 * np.pi
i_min = np.argmin(np.abs(w-1))
i_max = np.argmin(np.abs(w-4))

plt.figure()
plt.subplot(2,1,1)
plt.plot(w[i_min:i_max], np.abs(E_out[i_min:i_max]), 'r')
plt.ylabel('Spectral amplitude (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(w[i_min:i_max], phase_shifted[i_min:i_max], 'b')
plt.xlabel('Angular frequency (PHz)')
plt.ylabel('Spectral phase (rad)')
plt.grid()
plt.show()
plt.tight_layout()

intensity = np.abs(E_out) ** 2
GD = np.gradient(phase_shifted, dw)
GD0 = GD[indexw0]
relGD = GD-GD0

plt.figure()
plt.subplot(2,1,1)
plt.plot(w[i_min:i_max], np.abs(intensity[i_min:i_max]), 'r')
plt.ylabel('Spectral intensity (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(w[i_min:i_max], relGD[i_min:i_max], 'b')
plt.xlabel('Angular frequency (PHz)')
plt.ylabel('Rel. group delay (fs)')
plt.grid()
plt.show()
plt.tight_layout()

imax = np.argmax(intensity)
i0 = np.argmin(np.abs(intensity[:imax] - max(intensity)/2))
i1 = np.argmin(np.abs(intensity[imax:] - max(intensity)/2)) + imax
w1 = w[i0]
w2 = w[i1]
delta_GD = np.abs(GD[i0] - GD[i1])
print('Group delay difference (fused silica) = {0:0.3f} fs'.format(delta_GD))
print('GD(w0) (fused silica) = {0:0.3f} fs'.format(GD0))


GDD = np.gradient(GD, dw) # Group delay dispersion [fs^2]
GDD0 = GDD[indexw0]
print('GDD(w0) (fused silica) = {0:0.2f} fs^2'.format(GDD0))
# This value is used to estimate the pulse duration of the laser pulse
estimated_duration = tau_0 * np.sqrt(1 + (4 * np.log(2) * GDD0/(tau_0 ** 2))**2)
print('Duration from formula (fused silica) = {0:0.3f} fs'.format(estimated_duration))

Ew0 = np.fft.fftshift(E_out)
E_t_0 = np.fft.ifft(Ew0) * (omega / (2 * np.pi))
t0 = np.fft.fftfreq(len(Ew0), dw / (2 * np.pi))
Et = np.fft.fftshift(E_t_0)
t = np.fft.fftshift(t0)

fi_t = np.unwrap(np.angle(Et))
f0 = np.argmin(np.abs(t - GD0))
n = int(fi_t[f0]/(2 * np.pi))
fi_correct = fi_t - n * 2 * np.pi

t_min = np.argmin(np.abs(t - GD0 + 2 * estimated_duration))
t_max = np.argmin(np.abs(t - GD0 - 2 * estimated_duration))

plt.figure()
plt.subplot(2,1,1)
plt.plot(t[t_min:t_max], np.abs(Et[t_min:t_max]), 'r--')
plt.plot(t[t_min:t_max], np.real(Et[t_min:t_max]), 'r')
plt.ylabel('Electric field (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(t[t_min:t_max], fi_correct[t_min:t_max], 'b')
plt.xlabel('Time (fs)')
plt.ylabel('Temporal phase (rad)')
plt.grid()
plt.show()
plt.tight_layout()

intensity_t = np.abs(Et) ** 2
dt = 2 * np.pi / omega
w_inst_t = np.gradient(fi_correct, dt)

plt.figure()
plt.subplot(2,1,1)
plt.plot(t[t_min:t_max], intensity_t[t_min:t_max], 'r')
plt.ylabel('Intensity (a.u.)')
plt.title('I(t) and w_ins(t) of the output pulse (Fused Silica)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(t[t_min:t_max], w_inst_t[t_min:t_max], 'b')
plt.xlabel('Time (fs)')
plt.ylabel('Inst. ang. freq. (PHz)')
plt.grid()
plt.show()
plt.tight_layout()

itmax = np.argmax(np.abs(intensity_t))
it0 = np.argmin(np.abs(intensity_t[:itmax] - max(intensity_t)/2))
it1 = np.argmin(np.abs(intensity_t[itmax:] - max(intensity_t)/2)) + itmax
t1 = t[it0]
t2 = t[it1]
pulse_duration = t2 - t1
print("Duration from I(t) (Fused Silica) = {0:0.3f} fs".format(pulse_duration))


# Same exercises for BK7
n_w_BK7 = [refrindex(i, 'BK7') for i in w]
E_out = np.exp(-1j * (w/c) * n_w_BK7 * L) * Ew

spectral_phase = np.unwrap(-np.angle(E_out))
indexw0 = np.argmin(np.abs(w-w0))
n = int(spectral_phase[indexw0]/(2 * np.pi))
phase_shifted = spectral_phase  - n * 2 * np.pi

i_min = np.argmin(np.abs(w - 1))
i_max = np.argmin(np.abs(w - 4))

plt.figure()
plt.subplot(2,1,1)
plt.plot(w[i_min:i_max], np.abs(E_out[i_min:i_max]), 'r')
plt.ylabel('Spectral amplitude (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(w[i_min:i_max], phase_shifted[i_min:i_max], 'b')
plt.xlabel('Angular frequency (PHz)')
plt.ylabel('Spectral phase (rad)')
plt.grid()
plt.show()
plt.tight_layout()


intensity = np.abs(E_out) ** 2 
GD = np.gradient(phase_shifted, dw)
GD0 = GD[indexw0]
relGD = GD - GD0 

plt.figure()
plt.subplot(2,1,1)
plt.plot(w[i_min:i_max], np.abs(intensity[i_min:i_max]), 'r')
plt.ylabel('Spectral intensity (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(w[i_min:i_max], relGD[i_min:i_max], 'b')
plt.xlabel('Angular frequency (PHz)')
plt.ylabel('Rel. group delay (fs)')
plt.grid()
plt.show()
plt.tight_layout()

imax = np.argmax(intensity)
i0 = np.argmin(np.abs(intensity[:imax] - max(intensity)/2))
i1 = np.argmin(np.abs(intensity[imax:] - max(intensity)/2)) + imax
w1 = w[i0]
w2 = w[i1]
print('i0 = ', i0)
print('i1 = ', i1)
delta_GD = np.abs(GD[i0] - GD[i1])
print('Group delay difference (BK7) = {0:0.3f} fs'.format(delta_GD))
print('GD(w0) (BK7) = {0:0.3f} fs'.format(GD0))

# ---------------------------------------
# Exercise: Estimated pulse duration.
GDD = np.gradient(GD, dw) # Group delay dispersion [fs^2]
GDD0 = GDD[indexw0]
print('GDD(w0) (BK7) = {0:0.2f} fs^2'.format(GDD0))
estimated_duration = tau_0 * np.sqrt(1 + (4 * np.log(2) * GDD0/(tau_0 ** 2))**2)
print('Duration from formula (BK7) = {0:0.3f} fs'.format(estimated_duration))

# --------------------------------------------------------
# Exercise: E(t) and phi(t) of the output pulse
Ew0 = np.fft.fftshift(E_out)
E_t_0 = np.fft.ifft(Ew0) * (omega / (2 * np.pi))
t0 = np.fft.fftfreq(len(Ew0), dw / (2 * np.pi))
Et = np.fft.fftshift(E_t_0)
t = np.fft.fftshift(t0)

fi_t = np.unwrap(np.angle(Et))
f0 = np.argmin(np.abs(t - GD0))
n = int(fi_t[f0]/(2 * np.pi))
fi_correct = fi_t - n * 2 * np.pi

t_min = np.argmin(np.abs(t - GD0 + 2.5 * estimated_duration))
t_max = np.argmin(np.abs(t - GD0 - 2.5 * estimated_duration))

plt.figure()
plt.subplot(2,1,1)
plt.plot(t[t_min:t_max], np.abs(Et[t_min:t_max]), 'r--')
plt.plot(t[t_min:t_max], np.real(Et[t_min:t_max]), 'r')
plt.ylabel('Electric field (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(t[t_min:t_max], fi_correct[t_min:t_max], 'b')
plt.xlabel('Time (fs)')
plt.ylabel('Temporal phase (rad)')
plt.grid()
plt.show()
plt.tight_layout()

intensity_t = np.abs(Et) ** 2
dt = 2 * np.pi / omega
w_inst_t = np.gradient(fi_correct, dt)

plt.figure()
plt.subplot(2,1,1)
plt.plot(t[t_min:t_max], intensity_t[t_min:t_max], 'r')
plt.ylabel('Intensity (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(t[t_min:t_max], w_inst_t[t_min:t_max], 'b')
plt.xlabel('Time (fs)')
plt.ylabel('Inst. ang. freq. (PHz)')
plt.grid()
plt.show()
plt.tight_layout()

itmax = np.argmax(np.abs(intensity_t))
it0 = np.argmin(np.abs(intensity_t[:itmax] - max(intensity_t)/2))
it1 = np.argmin(np.abs(intensity_t[itmax:] - max(intensity_t)/2)) + itmax
t1 = t[it0]
t2 = t[it1]
pulse_duration = t2 - t1
print("Duration from I(t) (BK7) = {0:0.3f} fs".format(pulse_duration))


# Same exercises for SF10
n_w_SF10 = [refrindex(i, 'SF10') for i in w]
E_out = np.exp(-1j * (w/c) * n_w_SF10 * L) * Ew

spectral_phase = np.unwrap(-np.angle(E_out))
indexw0 = np.argmin(np.abs(w-w0))
n = int(spectral_phase[indexw0]/(2 * np.pi))
phase_shifted = spectral_phase  - n * 2 * np.pi

i_min = np.argmin(np.abs(w - 1))
i_max = np.argmin(np.abs(w - 4))

plt.figure()
plt.subplot(2,1,1)
plt.plot(w[i_min:i_max], np.abs(E_out[i_min:i_max]), 'r')
plt.ylabel('Spectral amplitude (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(w[i_min:i_max], phase_shifted[i_min:i_max], 'b')
plt.xlabel('Angular frequency (PHz)')
plt.ylabel('Spectral phase (rad)')
plt.grid()
plt.show()
plt.tight_layout()

intensity = np.abs(E_out) ** 2
GD = np.gradient(phase_shifted, dw)
GD0 = GD[indexw0]
relGD = GD - GD0


plt.figure()
plt.subplot(2,1,1)
plt.plot(w[i_min:i_max], np.abs(intensity[i_min:i_max]), 'r')
plt.ylabel('Spectral intensity (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(w[i_min:i_max], relGD[i_min:i_max], 'b')
plt.xlabel('Angular frequency (PHz)')
plt.ylabel('Rel. group delay (fs)')
plt.grid()
plt.show()
plt.tight_layout()

imax = np.argmax(intensity)
i0 = np.argmin(np.abs(intensity[:imax] - max(intensity)/2))
i1 = np.argmin(np.abs(intensity[imax:] - max(intensity)/2)) + imax
w1 = w[i0]
w2 = w[i1]
print('i0 = ', i0)
print('i1 = ', i1)
delta_GD = np.abs(GD[i0] - GD[i1])
print('Group delay difference (SF10) = {0:0.3f} fs'.format(delta_GD))
print('GD(w0) (SF10) = {0:0.3f} fs'.format(GD0))

GDD = np.gradient(GD, dw)
GDD0 = GDD[indexw0]
print('GDD(w0) (SF10) = {0:0.2f} fs^2'.format(GDD0))
estimated_duration = tau_0 * np.sqrt(1 + (4 * np.log(2) * GDD0/(tau_0 ** 2))**2)
print('Duration from formula (SF10) = {0:0.3f} fs'.format(estimated_duration))

Ew0 = np.fft.fftshift(E_out)
E_t_0 = np.fft.ifft(Ew0) * (omega / (2 * np.pi))
t0 = np.fft.fftfreq(len(Ew0), dw / (2 * np.pi))
Et = np.fft.fftshift(E_t_0)
t = np.fft.fftshift(t0)

fi_t = np.unwrap(np.angle(Et))
f0 = np.argmin(np.abs(t - GD0))
n = int(fi_t[f0]/(2 * np.pi))
fi_correct = fi_t - n * 2 * np.pi

t_min = np.argmin(np.abs(t - GD0 + 2.5 * estimated_duration))
t_max = np.argmin(np.abs(t - GD0 - 2.5 * estimated_duration))

plt.figure()
plt.subplot(2,1,1)
plt.plot(t[t_min:t_max], np.abs(Et[t_min:t_max]), 'r--')
plt.plot(t[t_min:t_max], np.real(Et[t_min:t_max]), 'r')
plt.ylabel('Electric field (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(t[t_min:t_max], fi_correct[t_min:t_max], 'b')
plt.xlabel('Time (fs)')
plt.ylabel('Temporal phase (rad)')
plt.grid()
plt.show()
plt.tight_layout()

intensity_t = np.abs(Et) ** 2
dt = 2 * np.pi / omega
w_inst_t = np.gradient(fi_correct, dt)


plt.figure()
plt.subplot(2,1,1)
plt.plot(t[t_min:t_max], intensity_t[t_min:t_max], 'r')
plt.ylabel('Intensity (a.u.)')
plt.grid()
plt.show()
plt.subplot(2,1,2)
plt.plot(t[t_min:t_max], w_inst_t[t_min:t_max], 'b')
plt.xlabel('Time (fs)')
plt.ylabel('Inst. ang. freq. (PHz)')
plt.grid()
plt.show()
plt.tight_layout()

itmax = np.argmax(np.abs(intensity_t))
it0 = np.argmin(np.abs(intensity_t[:itmax] - max(intensity_t)/2))
it1 = np.argmin(np.abs(intensity_t[itmax:] - max(intensity_t)/2)) + itmax
t1 = t[it0]
t2 = t[it1]
pulse_duration = t2 - t1
print("Duration from I(t) (SF10) = {0:0.3f} fs".format(pulse_duration))


