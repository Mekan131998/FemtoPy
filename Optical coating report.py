# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 07:41:02 2023

@author: Owner
"""

import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 

plt.close('all')

#49.83
ellip_BK7_49=pd.read_excel('Data/Optical coating/BK7_ellipsometry.xlsx', sheet_name="deg49.83", header=None, names=["E_eV", "Tan_Psi", "Cos_delta", "Psi", "Delta", "De_pol"])[2:]
ellip_Nb=pd.read_excel('Data/Optical coating/Nb2O5_fotometer.xlsx', header=None, names=["E_eV", "Tan_Psi", "Cos_delta", "Psi", "Delta", "De_pol"])[2:]

data_49=ellip_BK7_49[:-1].astype('float')

#54.83
ellip_BK7_54=pd.read_excel('Data/Optical coating/BK7_ellipsometry.xlsx', sheet_name="deg54.83", header=None, names=["E_eV", "Tan_Psi", "Cos_delta", "Psi", "Delta", "De_pol"])[2:]
data_54=ellip_BK7_54[:-1].astype('float')

#59.83
ellip_BK7_59=pd.read_excel('Data/Optical coating/BK7_ellipsometry.xlsx', sheet_name="deg59.83", header=None, names=["E_eV", "Tan_Psi", "Cos_delta", "Psi", "Delta", "De_pol"])[2:]
data_59=ellip_BK7_59[:-1].astype('float')



ref_BK7=pd.read_csv('Data/Optical coating/BK7-R.Sample.Raw.csv', header=None, names=["nm", "R"])[1:].astype('float')
tr_BK7=pd.read_csv('Data/Optical coating/BK7.Sample.Raw.csv', names=["nm", "T"])[1:].astype('float')
ref_FS=pd.read_csv('Data/Optical coating/FS-R.Sample.Raw.csv', names=["nm", "R"])[1:].astype('float')
tr_FS=pd.read_csv('Data/Optical coating/FS.Sample.Raw.csv', names=["nm", "T"])[1:].astype('float')

Tr=pd.read_excel('Data/Optical coating/Nb2O5_fotometer.xlsx', names=["nm", "T", "upper", "lower"]).astype('float')

# =============================================================================
# Plotting transmittance and reflectance on the same graph
# =============================================================================
fig, ax2 = plt.subplots()

color = 'tab:blue'

ax2.set_ylabel('Reflectance %', color=color) # we already handled the x-label with ax1
ax2.set_xlabel('Wavelength (nm)')
ax2.plot(ref_BK7["nm"], ref_BK7["R"], color='blue', label="BK7 ref")
ax2.plot(ref_FS["nm"], ref_FS["R"], color='green', label="FS ref")
ax2.tick_params(axis='y', labelcolor=color)
ax2.grid()
ax2.legend(loc='lower right')

ax1 = ax2.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'

ax1.plot(tr_BK7["nm"], tr_BK7["T"], color='red', label="BK7 tr")
ax1.plot(tr_FS["nm"], tr_FS['T'], color='purple', label="FS tr")
ax1.set_ylabel('Transmittance %', color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.legend(loc='center right')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('images\optical_coating_images/Reflectance_Transmittance.png')
plt.show()

l_BK7=6e6      #nm
l_FS=1e6       #nm

ref_BK7["R"]=ref_BK7["R"]*1e-2    #   % to dimensionless
ref_FS["R"]=ref_FS["R"]*1e-2      #   % to dimensionless 


k=-(ref_BK7["nm"]/(4*np.pi*l_BK7))*np.log(tr_BK7["T"]*1e-2/(1-ref_BK7["R"])**2)

k2=-(ref_FS["nm"]/(4*np.pi*l_FS))*np.log(tr_FS["T"]*1e-2/(1-ref_FS["R"])**2)


# n_tBK7=(np.sqrt(ref_BK7["R"])+1)/(1-np.sqrt(ref_BK7["R"]))
# n_tFS=(np.sqrt(ref_FS["R"])+1)/(1-np.sqrt(ref_FS["R"]))
r=(ref_BK7["R"]+1)/(ref_BK7["R"]-1)
n_tBK7=-r+np.sqrt(r**2-(k**2+1))

r1=(ref_FS["R"]+1)/(ref_FS["R"]-1)
n_tFS=-r1+np.sqrt(r1**2-(k2**2+1))

def n_FS(x0):
    x=x0*1e-3
    return (1+0.6961663/(1-(0.0684043/x)**2)+0.4079426/(1-(0.1162414/x)**2)+0.8974794/(1-(9.896161/x)**2))**.5

def n_BK7(x0):
    x=x0*1e-3
    return (1+1.03961212/(1-0.00600069867/x**2)+0.231792344/(1-0.0200179144/x**2)+1.01046945/(1-103.560653/x**2))**.5

n_BK_theo=n_BK7(ref_BK7["nm"])
n_FS_theo=n_FS(ref_FS["nm"])
# =============================================================================
#BK7 Plotting refractive index and the extinction coefficient on the same graph 
# =============================================================================
fig, ax2 = plt.subplots()

color = 'tab:blue'

ax2.set_ylabel('refractive index', color=color)  # we already handled the x-label with ax1
ax2.plot(ref_BK7["nm"], n_tBK7, color="blue", label='n of BK7')
ax2.plot(ref_BK7["nm"], n_BK_theo, color="black", linestyle="dashed", label='literature ', alpha=0.5)
ax2.tick_params(axis='y', labelcolor=color)
ax2.grid()
ax2.legend()

ax1 = ax2.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'

ax1.set_xlabel('wavelength, nm')
ax1.set_ylabel("Extinction coefficient", color=color)
ax1.plot(ref_BK7["nm"], k, color=color, label='k of BK7')
ax1.tick_params(axis='y', labelcolor=color)
ax1.legend(loc='upper center')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('images\optical_coating_images/k_ref_BK7.png')
plt.show()

plt.figure()
plt.plot(ref_FS["nm"], k2)
plt.show()

# =============================================================================
# FS 
# =============================================================================
fig, ax2 = plt.subplots()

color = 'tab:blue'

ax2.set_ylabel('refractive index', color=color)  # we already handled the x-label with ax1
ax2.plot(ref_FS["nm"], n_tFS, color="blue", label='n of FS')
ax2.plot(ref_FS["nm"], n_FS_theo, color="black", linestyle="dashed", label='literature ', alpha=0.5)
ax2.tick_params(axis='y', labelcolor=color)
ax2.grid()
ax2.legend()

ax1 = ax2.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'

ax1.set_xlabel('wavelength, nm')
ax1.set_ylabel("Extinction coefficient", color=color)
ax1.plot(ref_FS["nm"], k2, color=color, label='k of FS')
ax1.tick_params(axis='y', labelcolor=color)
ax1.legend(loc='upper center')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('images\optical_coating_images/k_ref_FS.png')
plt.show()

plt.figure()
plt.plot(ref_FS["nm"], k2)
plt.show()

# =============================================================================
# Transmission curve: 
# =============================================================================

T_L= Tr["lower"]*1e-2

#673.8 nm 
#550 nm 

plt.figure()
plt.plot(Tr["nm"], Tr["T"]*1e-2)
plt.plot(Tr["nm"], Tr["upper"]*1e-2, "g--", label="upper $T_U$")
plt.plot(Tr["nm"], T_L, "b--", label="lower $T_L$")
plt.vlines(470, 0, 0.885 , linestyles="dashed", color='grey', alpha=0.7)
plt.vlines(550, 0, 0.86 , linestyles="dashed", color='grey', alpha=0.7)
plt.annotate('470', xy=(470, 0.02), xytext=(470, 0.01), color='grey')
plt.annotate('550', xy=(552, 0.02), xytext=(552, 0.01), color='grey')
plt.xlabel("wavelength, nm ")
plt.ylabel("Transmission")
plt.grid()
plt.legend()
plt.savefig('images\optical_coating_images/Trans_layer_Upper_Lower.png')
plt.show()




# =============================================================================
# Single layer characterization 
# =============================================================================
lmbda=ref_BK7["nm"]
indexes_in_range = lmbda[(lmbda >= 350) & (lmbda <= 650)].index

T_L_actual=T_L[indexes_in_range]
n_tBK7_actual=n_tBK7[indexes_in_range]

N=2*n_tBK7_actual/T_L_actual-(n_tBK7_actual**2+1)/2  

n_l=np.sqrt(N+np.sqrt(N**2-n_tBK7_actual**2))

plt.figure()
plt.plot(lmbda[indexes_in_range], n_l)
plt.ylabel("refractive index of the layer")
plt.xlabel("wavelength (nm)")
plt.vlines(470, 0, 2.478, linestyles="dashed", color='grey', alpha=0.7)
plt.vlines(550, 0, 2.41, linestyles="dashed", color='grey', alpha=0.7)
plt.hlines(2.478, 335, 470, linestyles='dashed', color='grey', alpha=0.7)
plt.hlines(2.41, 335, 550, linestyles='dashed', color='grey', alpha=0.7)
plt.annotate('470', xy=(470, 2.01 ), xytext=(470, 2.01), color='grey')
plt.annotate('550', xy=(552, 2.01), xytext=(552, 2.01), color='grey')
plt.grid()
plt.ylim(2 , 4)
plt.savefig('images\optical_coating_images/Ref_layer.png')
plt.show()

# =============================================================================
# Ellipsometry 
# =============================================================================
h=4.135
c=299.792 
        # 49.83, 54.83 and 59.83

def last(data, angle):
    waveln=h*c/data["E_eV"]
    
    inx_range=waveln[(waveln >= 250) & (waveln <= 950)].index
    Psi_degree=data["Psi"]
    Delta_degree=data["Delta"]
    Delta_rad=Delta_degree*np.pi/180
    n_lit=n_BK7(waveln[inx_range])
    
    plt.figure()
    plt.subplot(211)
    plt.title("Angle ${}^0$".format(angle))
    plt.plot(waveln[inx_range], Psi_degree[inx_range])
    plt.ylabel("$\Psi$")
    plt.xlabel("wavelength (nm)")
    plt.grid()
    
    plt.subplot(212)
    plt.plot(waveln[inx_range], Delta_degree[inx_range])
    plt.ylabel("$\Delta$")
    plt.xlabel("wavelength (nm)")
    plt.grid()
    plt.tight_layout()
    plt.savefig('images\optical_coating_images/psi_delta_{}.png'.format(angle))
    plt.show()
    
    phi=angle*np.pi/180
    
    ro=data["Tan_Psi"]*np.exp(1j*Delta_rad)
    
    e=(np.sin(phi))**2*(1+np.tan(phi)**2*((1-ro)/(1+ro))**2)
    
    n_complex=np.sqrt(e)
    n_ref=np.real(n_complex)
    k_ex=np.imag(n_complex)

    plt.figure()
    plt.subplot(211)
    plt.title("Angle ${}^0$".format(angle))
    plt.plot(waveln[inx_range], n_ref[inx_range], label='calculated')
    plt.plot(waveln[inx_range], n_lit, color="red", linestyle="dashed", label='literature')
    plt.ylabel("n")
    plt.xlabel("wavelength (nm)")
    plt.legend()
    plt.grid()
    
    plt.subplot(212)
    plt.plot(waveln[inx_range], k_ex[inx_range], label="calculated k")
    plt.hlines(9.7525e-9, min(waveln[inx_range]), max(waveln[inx_range]), linestyles="dashed", color="black", label="literature")
    plt.ylabel("k")
    plt.xlabel("wavelength (nm)")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('images\optical_coating_images/n_k_final_{}.png'.format(angle))
    plt.show()


# plotting for 49.83:
last(data_49, 49.83)

# plotting for 54.83
last(data_54, 54.83)

# plotting for 59.83
last(data_59, 59.83)

d=(470*550)/(2*(550*2.478-470*2.413))
print(d)


