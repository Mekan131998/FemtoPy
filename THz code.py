# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 18:08:05 2023

@author: Owner
"""

import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.integrate import cumulative_trapezoid

# =============================================================================
# 1. Pyrocamera data
# =============================================================================
plt.close('all')

inten_prfl=pd.read_csv('data/THz_spectro/Pyrocam_intensity_profile.csv', header=None)

Ix=np.sum(inten_prfl, axis=0)
Iy=np.sum(inten_prfl, axis=1)

Ix_n=Ix/max(Ix)
Iy_n=Iy/max(Iy)

# Create a 2x2 grid for subplots using GridSpec
fig = plt.figure(figsize=(8, 8))
gs = GridSpec(2, 2, width_ratios=[1, 4], height_ratios=[4, 1])

# Plot sum along y axis
ax_sum_y = fig.add_subplot(gs[0, 0], sharey=None)
ax_sum_y.plot(Iy_n, np.arange(len(Iy)))
ax_sum_y.set_xlabel('Sum along X axis')

# Plot 2D Gaussian laser beam intensity
ax_intensity = fig.add_subplot(gs[0, 1], sharey=ax_sum_y)
cmap = plt.get_cmap('viridis')
im = ax_intensity.imshow(inten_prfl, cmap=cmap)
# fig.colorbar(im, ax=ax_intensity, label='Intensity')

# Plot sum along x axis
ax_sum_x = fig.add_subplot(gs[1, 1], sharex=ax_intensity)
ax_sum_x.plot(Ix_n)
ax_sum_x.set_ylabel('Sum along Y axis')

# Remove x and y axis labels for the shared axes
plt.setp(ax_sum_y.get_yticklabels(), visible=False)
plt.setp(ax_sum_x.get_xticklabels(), visible=False)
plt.savefig("Data/THz_spectro/Beam_profile_THz.png")
plt.show()



# Voltage signal: 
angle=np.array([100, 90, 80, 70, 60, 50, 40, 30])
S=np.array([8.35, 22.096, 46.83, 83.372, 123.7, 159.28, 184.95, 195.39])  #in mV

T_pmp=0.934
T_cb=0.722 
T=T_cb*T_pmp 
R=55.16 

fr=1000 #Hz

W=S/(T*R*fr) #mV/Hz
W=W*1e3      # 1e-3 mV/Hz 
print("THz pulse energy values 1e-3 mV: {}".format(W))



# finding beam size: 

x=np.arange(len(Ix_n))
y=np.arange(len(Iy_n))

x1=x[Ix_n>max(Ix_n/2)/2][0]
x2=x[Ix_n>max(Ix_n/2)/2][-1]

y1=y[Iy_n>max(Iy_n/2)/2][0]
y2=y[Iy_n>max(Iy_n/2)/2][-1]

d1=x2-x1
diam1=1.7*d1

d2=y2-y1
diam2=1.7*d2
# ratio between pixel and size:
    
r=25.6/320      # in mm 

diam1=diam1*r
diam2=diam2*r

# Area:
A=np.pi*diam1*diam2/4

print("d1[mm]= {}".format(diam1))
print("d2[mm]= {}".format(diam2))
print("Area[mm2]= {0:0.2f}".format(A))



# =============================================================================
# 2. Waveform measurements 
# =============================================================================

air_30=pd.read_csv('data/THz_spectro/laserEnergy_3.58mJ_air_25ps_wg30_16um_st3_sam10_RT_delay_6.07.txt', sep="\t")
Ge_30=pd.read_csv('data/THz_spectro/laserEnergy_3.58mJ_nGe_20ps_wg30_16um_st1_sam10_RT_N2_delay_6.07.txt', sep="\t")
Ge_90=pd.read_csv('data/THz_spectro/laserEnergy_3.58mJ_nGe_20ps_wg90_16um_st3_sam10_RT_N2_delay_6.07.txt', sep="\t")
ref_30=pd.read_csv('data/THz_spectro/laserEnergy_3.58mJ_REF_20ps_wg30_16um_st1_sam10_RT_N2_delay_6.07.txt', sep="\t")
ref_90=pd.read_csv('data/THz_spectro/laserEnergy_3.58mJ_REF_20ps_wg90_16um_st3_sam10_RT_N2_delay_6.07.txt', sep="\t")


# calculating integral: 
air_30_sq=air_30["Norm[arb.u.]"]**2 
Ge_30_sq=Ge_30["Norm[arb.u.]"]**2 
Ge_90_sq=Ge_90["Norm[arb.u.]"]**2 
ref_30_sq=ref_30["Norm[arb.u.]"]**2 
ref_90_sq=ref_90["Norm[arb.u.]"]**2 

i_ref_30=cumulative_trapezoid(ref_30_sq, ref_30["Time[ps]"], initial=0)
i_ref_90=cumulative_trapezoid(ref_90_sq, ref_90["Time[ps]"], initial=0)
i_Ge_30=cumulative_trapezoid(Ge_30_sq, Ge_30["Time[ps]"], initial=0)
i_Ge_90=cumulative_trapezoid(Ge_90_sq, Ge_90["Time[ps]"], initial=0)


print("Integral for ref_30 : {0:0.4f} \n".format(max(i_ref_30)))
print("Integral for ref_90 : {0:0.4f} \n".format(max(i_ref_90)))
print("Integral for Ge_30 : {0:0.4f} \n".format(max(i_Ge_30)))
print("Integral for Ge_90 : {0:0.4f} \n".format(max(i_Ge_90)))

integ_ref_30=max(i_ref_30)
integ_ref_90=max(i_ref_90)
integ_Ge_30=max(i_Ge_30)
integ_Ge_90=max(i_Ge_90)
              
# Energy calculation: SI units:
e0=8.85*1e-12
c=3*1e8 

#in 1e-3 mV dimension:

W30=5.2528
W90=0.594

E0_ref_30=np.sqrt(W30/(e0*c*A*integ_ref_30))
E0_ref_90=np.sqrt(W90/(e0*c*A*integ_ref_90))
E0_Ge_90=np.sqrt(W90/(e0*c*A*integ_Ge_90))
E0_Ge_30=np.sqrt(W30/(e0*c*A*integ_Ge_30))

print("E0_ref_30 [V/m] ={}\n".format(E0_ref_30))
print("E0_ref_90 [V/m] ={}\n".format(E0_ref_90))
print("E0_Ge_90 [V/m] ={}\n".format(E0_Ge_90))
print("E0_Ge_30 [V/m] ={}\n".format(E0_Ge_30))

C_ref_30=E0_ref_30*ref_30["AVG[arb.u.]"]
C_ref_90=E0_ref_90*ref_90["AVG[arb.u.]"]
C_Ge_30=E0_Ge_30*Ge_30["AVG[arb.u.]"]
C_Ge_90=E0_Ge_90*Ge_90["AVG[arb.u.]"]

# Norm[arb.u.]
# Time[ps]

plt.figure()
plt.plot(Ge_30["Time[ps]"],C_Ge_30, color='blue', label="sample pulse Ge @30")
plt.plot(ref_30["Time[ps]"], C_ref_30, color="green", label="reference pulse @30")
plt.xlabel("Time delay, (ps)")
plt.ylabel("THz field strength (a.u)")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('Data/THz_spectro/norm_Ge30_vs_ref30.png')
plt.show()

plt.figure()
plt.plot(Ge_90["Time[ps]"], C_Ge_90, color='blue', label="sample pulse Ge @90")
plt.plot(ref_90["Time[ps]"], C_ref_90, color="green", label="reference pulse @ 90")
plt.xlabel("Time delay, (ps)")
plt.ylabel("THz field strength (a.u)")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('Data/THz_spectro/Ge90_vs_ref90.png')
plt.show()

# =============================================================================
# Fourier transform 
# =============================================================================

def fourierTransform(Et, k_fac):
    T=max(Et["Time[ps]"])-min(Et["Time[ps]"])
    N=2*16
    dt=T/N
    Et_0T=np.fft.fftshift(k_fac*Et["AVG[arb.u.]"])
    C_0OM=np.fft.fft(Et_0T)/N
    freq_0OM=np.fft.fftfreq(len(Et_0T), dt)
    C=np.fft.fftshift(C_0OM)
    freq=np.fft.fftshift(freq_0OM)
    om=2*np.pi*freq
    return C, om

ft_Ge30, om_G30=fourierTransform(Ge_30, E0_Ge_30)
ft_ref30, om_r30=fourierTransform(ref_30, E0_ref_30)




ft_Ge90, om_G90=fourierTransform(Ge_90, E0_Ge_90)
ft_ref90, om_r90=fourierTransform(ref_90, E0_ref_90)



T_w_30=ft_Ge30/ft_ref30
T_w_90=ft_Ge90/ft_ref90

d=0.5   #in mm 

ref_index_30=1+np.angle(T_w_30)*c/(om_G30*d*1e12)
alpha_30=-(2/d)*np.log((ref_index_30+1)**2/(4*ref_index_30)*abs(T_w_30))

ref_index_90=1+np.angle(T_w_90)*c/(om_G90*d*1e12)
alpha_90=-(2/d)*np.log((ref_index_90+1)**2/(4*ref_index_90)*abs(T_w_30))

plt.figure(figsize=(6,18))
plt.subplot(311)
plt.title("a)")
plt.plot(om_G30, abs(ft_Ge30), color="blue", label="Ge @30")
plt.plot(om_r30, abs(ft_ref30), color='green', label="ref @30")
plt.ylabel("THz field strength, (a.u)")
plt.xlim(0, max(om_G30))
plt.yscale("log")
plt.grid()
plt.legend()
# plt.savefig('Data/THz_spectro/freq_Ge30_vs_ref30.png')
# plt.show()

plt.subplot(312)
plt.title('b)')
plt.plot(om_G30, alpha_30, color="blue", label="Absorbance @30")
plt.ylabel("Absorbance, (a.u)")
plt.xlim(0, max(om_G30))
plt.grid()
plt.legend()
# plt.savefig('Data/THz_spectro/Absorbance_30.png')
# plt.show()

plt.subplot(313)
plt.title("c)")
plt.plot(om_G30, ref_index_30, color="blue", label="ref. index @30")
plt.xlabel("angular frequency, (THz)")
plt.ylabel("Refractive index, (a.u)")
plt.xlim(0, max(om_G30))
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('Data/THz_spectro/Ref_index_30.png')
plt.show()



plt.figure(figsize=(6,18))
plt.subplot(311)
plt.title("a)")
plt.plot(om_G90, abs(ft_Ge90), color="blue", label="Ge @90")
plt.plot(om_r90, abs(ft_ref90), color='green', label="ref @90")
# plt.xlabel("angular frequency, (THz)")
plt.ylabel("THz field strength, (a.u)")
plt.xlim(0, max(om_G90))
plt.yscale("log")
plt.grid()
plt.legend()
# plt.savefig('Data/THz_spectro/freq_Ge90_vs_ref90.png')
# plt.show()

plt.subplot(312)
plt.title("b)")
plt.plot(om_G90, alpha_90, color="blue", label="absorbance @90")
# plt.xlabel("angular frequency, (THz)")
plt.ylabel("Absorbance, (a.u)")
plt.xlim(0, max(om_G90))
plt.grid()
plt.legend()
# plt.savefig('Data/THz_spectro/Absorbance_90.png')


plt.subplot(313)
plt.title("c)")
plt.plot(om_G90, ref_index_90, color="blue", label="ref. index @90")
plt.xlabel("angular frequency, (THz)")
plt.ylabel("Refractive index, (a.u)")
plt.xlim(0, max(om_G90))
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('Data/THz_spectro/Ref_index_90.png')
plt.show()

