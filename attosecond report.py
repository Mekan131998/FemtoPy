# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 18:19:22 2023

@author: Owner
"""

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from matplotlib.gridspec import GridSpec
import pandas as pd 
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.close('all')

image_path = 'Data\Applications_atto/Task11_holey_mirror.png'
img = Image.open(image_path)
img.load()
data = np.asarray(img, dtype="int32")

sum_x=np.sum(data, axis=0)
sum_y=np.sum(data, axis=1)

norm_x=sum_x/max(sum_x)
norm_y=sum_y/max(sum_y)

len_x=np.arange(len(sum_x))
len_y=np.arange(len(sum_y))

# handing outliers: 1244 1488

def outliers(array):
    mod_array=array.copy()
    mod_array[1244:1488]=np.nan       
    return mod_array

def outliers2(array):
    mod_array=array.copy()
    mod_array[531:737]=np.nan       
    return mod_array

norm_x_clean=outliers(norm_x)
norm_y_clean=outliers2(norm_y)
ratio=7/165 

# curve fitting: 
# define the model function: 
def model_f(x, w_in, x0, A, B):
    return A*np.exp(-2*(x-x0)**2/(w_in**2))+B 

valid_indices = ~np.isnan(norm_x_clean)

# curve fit: 
popt, pcov=curve_fit(model_f, len_x[valid_indices], norm_x[valid_indices], p0=[221, 1378, 1, 0])
w_opt, x0_opt, A_opt, B_opt=popt

x_model=np.arange(min(len_x), max(len_x)+1, 1)
y_model=model_f(x_model, w_opt, x0_opt, A_opt, B_opt)

valid_indices2 = ~np.isnan(norm_y_clean)

popt2, pcov2=curve_fit(model_f, len_y[valid_indices2], norm_y[valid_indices2], p0=[221, 500, 1, 0])
w_opt2, x0_opt2, A_opt2, B_opt2=popt2



x_model2=np.arange(min(len_y), max(len_y)+1, 1)
y_model2=model_f(x_model2, w_opt2, x0_opt2, A_opt2, B_opt2)

# =============================================================================
# plt.figure()
# plt.imshow(data)
# plt.show()
# =============================================================================

# =============================================================================
# 
# plt.figure()
# plt.subplot(211)
# plt.scatter(len_x*ratio, norm_x_clean, edgecolors=None, s=3, c='blue', label="data")
# plt.plot(x_model*ratio, y_model, "g--", label="fitting")
# plt.xlabel("distance mm")
# plt.ylabel("Intensity (arb. units)")
# plt.legend()    
# 
# plt.subplot(212)
# plt.scatter(len_y*ratio, norm_y_clean, edgecolors=None, s=3,  c='blue', label='data')
# plt.plot(x_model2*ratio, y_model2, "r--", label='fitting')
# plt.xlabel("distance (mm)")
# plt.ylabel("Intensity (arb. units)")
# plt.legend()
# plt.tight_layout()
# plt.show()
# =============================================================================

fig = plt.figure(figsize=(8, 8))
gs = GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[4, 1])

# Plot sum along y axis
# Plot 2D Gaussian laser beam intensity
ax_intensity = fig.add_subplot(gs[0, 0])
cmap = plt.get_cmap('viridis')
im = ax_intensity.imshow(data, cmap=cmap)
# fig.colorbar(data, ax=ax_intensity, label='Intensity')

ax_sum_y = fig.add_subplot(gs[0, 1], sharey=ax_intensity)
ax_sum_y.scatter(norm_y_clean, np.arange(len(norm_y_clean)), edgecolors=None, s=3, c='blue', label="data")
ax_sum_y.plot(y_model2, np.arange(len(y_model2)), linestyle="dashed", color="red", label="fitting")
ax_sum_y.set_xlabel('Sum along X axis')

# Plot sum along x axis
ax_sum_x = fig.add_subplot(gs[1, 0], sharex=ax_intensity)
ax_sum_x.scatter(np.arange(len(norm_x_clean)), norm_x_clean, edgecolors=None, s=3, c='blue', label="data")
ax_sum_x.plot(y_model, linestyle="dashed", label="fitting")
ax_sum_x.set_ylabel('Sum along Y axis')
ax_sum_x.set_xlabel('pixel')

# Remove x and y axis labels for the shared axes
# plt.setp(ax_intensity.get_xticklabels(), visible=False)
# plt.setp(ax_intensity.get_yticklabels(), visible=False)
plt.setp(ax_sum_y.get_yticklabels(), visible=False)
# plt.setp(ax_sum_x.get_xticklabels(), visible=False)

ax_sum_y.legend(loc='upper left')
ax_sum_x.legend(loc='upper left')
plt.savefig('images/images_atto/Gaussian_fitting.png')

plt.show()

# =============================================================================
# 12. Cutoff energy calculation
# =============================================================================
w_in=7.2    #mm 
lmbda0=1030   #in nm 
f=900         #in mm 
w=lmbda0*f/(np.pi*w_in)      #in nm 

w1=w*1e-7    #in sm 
tau=40e-15        #in s 
E_n=1e-3          # J 

I0=2*E_n*np.sqrt((4*np.log(2))/np.pi)/(tau*w1**2*np.pi)

lmbda=lmbda0*1e-3    #in um

r=w/2

I_L=I0*np.exp(-2*r**2/w**2)
I_L1=I0*(1-0.4)

Up=9.33*1e-14*I_L1*lmbda**2    # in eV 

Ip=15.76   #eV 

E_cut=Ip+3.2*Up

Up0=9.33*1e-14*I0*lmbda**2    # in eV

E_cut0=Ip+3.2*Up0



# =============================================================================
# 13
# =============================================================================

img0 = Image.open('Data\Applications_atto/Task13_nofilter.png')
img0.load()
no_filter_data = np.asarray(img0, dtype="int32")


img1 = Image.open('Data\Applications_atto/Task13_Mgfilter.png')
img1.load()
mg_filter_data = np.asarray(img1, dtype="int32")


img2 = Image.open('Data\Applications_atto/Task13_Tefilter.png')
img2.load()
te_filter_data = np.asarray(img2, dtype="int32")
# =============================================================================
# 
# plt.figure()
# plt.title("no filter")
# plt.imshow(no_filter_data)
# plt.xlabel("pixels")
# plt.ylabel("pixels")
# plt.show()
# 
# plt.figure()
# plt.title('Mg filter')
# plt.imshow(mg_filter_data)
# plt.xlabel("pixels")
# plt.ylabel("pixels")
# plt.show()
# 
# plt.figure()
# plt.title("Te filter")
# plt.imshow(te_filter_data)
# plt.xlabel("pixels")
# plt.ylabel("pixels")
# plt.show()
# 
# =============================================================================
#removing background: 

bgrd_no_filter=no_filter_data[0:396, 0:400]
bgrd_mg_filter=mg_filter_data[0:396, 0:400]
bgrd_te_filter=te_filter_data[0:396, 0:400]

# =============================================================================
# plt.figure()
# plt.subplot(131)
# plt.title("Bgrd no filter")
# plt.imshow(bgrd_no_filter)
# 
# plt.subplot(132)
# plt.title("Bgrd mg filter")
# plt.imshow(bgrd_mg_filter)
# 
# plt.subplot(133)
# plt.title("Bgrd te filter")
# plt.imshow(bgrd_te_filter)
# plt.show()
# =============================================================================

avg_no_filter=np.mean(np.mean(bgrd_no_filter, axis=1))
avg_mg_filter=np.mean(np.mean(bgrd_mg_filter, axis=1))
avg_te_filter=np.mean(np.mean(bgrd_te_filter, axis=1))

no_ftr_c=no_filter_data-avg_no_filter
mg_ftr_c=mg_filter_data-avg_mg_filter
te_ftr_c=te_filter_data-avg_te_filter

sumX_no=np.sum(no_ftr_c, axis=0)
sumY_no=np.sum(no_ftr_c, axis=1)

sumX_mg=np.sum(mg_ftr_c, axis=0)
sumY_mg=np.sum(mg_ftr_c, axis=1)

sumX_te=np.sum(te_ftr_c, axis=0)
sumY_te=np.sum(te_ftr_c, axis=1)

x=np.arange(1, len(sumX_no)+1, 1)
y=np.arange(1, len(sumY_no)+1, 1)
#converting pixels to mm 
x_mm=x*ratio 
y_mm=y*ratio 
# =============================================================================
# 
# plt.figure()
# plt.subplot(331)
# plt.title("no filter removed")
# plt.imshow(no_ftr_c)
# 
# plt.subplot(332)
# plt.title("mg filter removed")
# plt.imshow(mg_ftr_c)
# 
# plt.subplot(333)
# plt.title("te filter")
# plt.imshow(te_ftr_c)
# 
# plt.subplot(334)
# plt.plot(x, sumX_no)
# plt.xlabel("pixel")
# 
# plt.subplot(335)
# plt.plot(x, sumX_mg)
# plt.xlabel("pixel")
# 
# plt.subplot(336)
# plt.plot(x, sumX_te)
# plt.xlabel("pixel")
# 
# plt.subplot(337)
# plt.plot(y, sumY_no)
# plt.xlabel("pixel")
# 
# plt.subplot(338)
# plt.plot(y, sumY_mg)
# plt.xlabel("pixel")
# 
# plt.subplot(339)
# plt.plot(y, sumY_te)
# plt.xlabel("pixel")
# 
# plt.tight_layout()
# plt.show()
# =============================================================================

# =============================================================================
# Calibrated image
# =============================================================================

fig = plt.figure(figsize=(8, 8))
gs = GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[4, 1])

# Plot sum along y axis
dimx, dimy=np.shape(no_ftr_c)
ax_intensity = fig.add_subplot(gs[0, 0])
cmap = plt.get_cmap('viridis')
im = ax_intensity.imshow(no_ftr_c, cmap=cmap)
# fig.colorbar(data, ax=ax_intensity, label='Intensity')

# along y 
ax_sum_y = fig.add_subplot(gs[0, 1], sharey=ax_intensity)
ax_sum_y.plot(sumY_no, np.arange(len(sumY_no)), color="Green")
ax_sum_y.set_xlabel('Sum along X axis')

# Plot sum along x axis
ax_sum_x = fig.add_subplot(gs[1, 0], sharex=ax_intensity)
ax_sum_x.plot(sumX_no)
ax_sum_x.set_ylabel('Sum along Y axis')
ax_sum_x.set_xlabel('pixel')
ax_intensity.set_title('No filter FFS backgound subtracted image')
ax_intensity.hlines(y=187, xmin=622, xmax=1312, linewidth=2, color='r', linestyle="dashed")
ax_intensity.hlines(y=622, xmin=622, xmax=1312, linewidth=2, color='r', linestyle="dashed")
ax_intensity.vlines(x=622, ymin=187, ymax=622, linewidth=2, color='r', linestyle="dashed")
ax_intensity.vlines(x=1312, ymin=187, ymax=622, linewidth=2, color='r', linestyle="dashed")

# Remove x and y axis labels for the shared axes
# plt.setp(ax_intensity.get_xticklabels(), visible=False)
# plt.setp(ax_intensity.get_yticklabels(), visible=False)
plt.setp(ax_sum_y.get_yticklabels(), visible=False)
# plt.setp(ax_sum_x.get_xticklabels(), visible=False)

plt.savefig('images/images_atto/FFS_image_no_filter.png')

plt.show()

# =============================================================================
# for mg
# =============================================================================
fig = plt.figure(figsize=(8, 8))
gs = GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[4, 1])

# Plot sum along y axis
dimx, dimy=np.shape(no_ftr_c)
ax_intensity = fig.add_subplot(gs[0, 0])
cmap = plt.get_cmap('viridis')
im = ax_intensity.imshow(mg_ftr_c, cmap=cmap)
# fig.colorbar(data, ax=ax_intensity, label='Intensity')

# along y 
ax_sum_y = fig.add_subplot(gs[0, 1], sharey=ax_intensity)
ax_sum_y.plot(sumY_mg, np.arange(len(sumY_mg)), color="Blue")
ax_sum_y.set_xlabel('Sum along X axis')

# Plot sum along x axis
ax_sum_x = fig.add_subplot(gs[1, 0], sharex=ax_intensity)
ax_sum_x.plot(sumX_mg, color="orange")
ax_sum_x.set_ylabel('Sum along Y axis')
ax_sum_x.set_xlabel('pixel')
ax_intensity.set_title('Mg filter FFS backgound subtracted image')
ax_intensity.hlines(y=187, xmin=622, xmax=1312, linewidth=2, color='r', linestyle="dashed")
ax_intensity.hlines(y=622, xmin=622, xmax=1312, linewidth=2, color='r', linestyle="dashed")
ax_intensity.vlines(x=622, ymin=187, ymax=622, linewidth=2, color='r', linestyle="dashed")
ax_intensity.vlines(x=1312, ymin=187, ymax=622, linewidth=2, color='r', linestyle="dashed")

# Remove x and y axis labels for the shared axes
# plt.setp(ax_intensity.get_xticklabels(), visible=False)
# plt.setp(ax_intensity.get_yticklabels(), visible=False)
plt.setp(ax_sum_y.get_yticklabels(), visible=False)
# plt.setp(ax_sum_x.get_xticklabels(), visible=False)

plt.savefig('images/images_atto/FFS_image_mg.png')

plt.show()

# =============================================================================
# Te
# =============================================================================

fig = plt.figure(figsize=(8, 8))
gs = GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[4, 1])

# Plot sum along y axis
dimx, dimy=np.shape(no_ftr_c)
ax_intensity = fig.add_subplot(gs[0, 0])
cmap = plt.get_cmap('viridis')
im = ax_intensity.imshow(te_ftr_c, cmap=cmap)
# fig.colorbar(data, ax=ax_intensity, label='Intensity')

# along y 
ax_sum_y = fig.add_subplot(gs[0, 1], sharey=ax_intensity)
ax_sum_y.plot(sumY_te, np.arange(len(sumY_te)), color="Blue")
ax_sum_y.set_xlabel('Sum along X axis')

# Plot sum along x axis
ax_sum_x = fig.add_subplot(gs[1, 0], sharex=ax_intensity)
ax_sum_x.plot(sumX_te, color="green")
ax_sum_x.set_ylabel('Sum along Y axis')
ax_sum_x.set_xlabel('pixel')
ax_intensity.set_title('Te filter FFS backgound subtracted image')
ax_intensity.hlines(y=187, xmin=622, xmax=1312, linewidth=2, color='r', linestyle="dashed")
ax_intensity.hlines(y=622, xmin=622, xmax=1312, linewidth=2, color='r', linestyle="dashed")
ax_intensity.vlines(x=622, ymin=187, ymax=622, linewidth=2, color='r', linestyle="dashed")
ax_intensity.vlines(x=1312, ymin=187, ymax=622, linewidth=2, color='r', linestyle="dashed")
# Remove x and y axis labels for the shared axes
# plt.setp(ax_intensity.get_xticklabels(), visible=False)
# plt.setp(ax_intensity.get_yticklabels(), visible=False)
plt.setp(ax_sum_y.get_yticklabels(), visible=False)
# plt.setp(ax_sum_x.get_xticklabels(), visible=False)

plt.savefig('images/images_atto/FFS_image_te.png')

plt.show()


# =============================================================================
# Finding transmission T 
# =============================================================================
#spatially integrating the data:
#187:608 for no filter, for w: 622:1537
#314:608 for mg filter, for w: 724:1318
#186:622 for te filter, for w: 767:1312
# 622:1312


x_no=np.arange(622, 1312, 1)

#Spatially integrating the data: 
S_no=sumY_no[187:622]

S_mg=sumY_mg[187:622]

S_te=sumY_te[187:622]

#SpEctrally integrating the data: 

S_w_no=sumX_no[622:1312]

S_w_mg=sumX_mg[622:1312]

S_w_te=sumX_te[622:1312]


trans_mg=pd.read_csv("Data/Applications_atto/Mg_transmission.txt", sep="\t", names=["energy", "trans"])
# =============================================================================
# sumX_nos=np.sum(no_ftr_c[187:622, 622:1312], axis=0)
# sumX_Mgs=np.sum(mg_ftr_c[187:622, 622:1312], axis=0)
# 
# sumY_nos=np.sum(no_ftr_c[187:622, 622:1312], axis=1)
# sumY_Mgs=np.sum(mg_ftr_c[187:622, 622:1312], axis=1)
# 
# Ts=sumX_Mgs/sumX_nos
# 
# =============================================================================
T=S_mg/S_no

T_w=S_w_mg/S_w_no
T_te=S_w_te/S_w_no

plt.figure(figsize=(6, 8))
plt.subplot(211)
plt.plot(x[622:1312], sumX_no[622:1312]/2, label="no filter*")
plt.plot(x[622:1312], sumX_mg[622:1312], label="mg filter")
plt.plot(x[622:1312], sumX_te[622:1312], label="te filter")
plt.vlines(x=915, ymin=0, ymax=1e5, linestyle="dashed", color="red")
plt.xlabel("pixel")
plt.ylabel("Intensity (arb. units)")
plt.legend()
plt.grid()  

plt.subplot(212)
plt.plot(x_no, np.abs(T_w), label="Mg", color="blue")
plt.plot(x_no, np.abs(T_te), label="Te", color="orange")
plt.grid()
plt.xlabel("pixel")
plt.ylabel("Transmission")
plt.vlines(x=915, ymin=0, ymax=0.42, linestyle="dashed", color="red")
plt.legend()
plt.tight_layout()
plt.savefig('images/images_atto/measured_Mg_Te_trans.png')
plt.show()

trans_te=pd.read_csv("Data/Applications_atto/Te_transmission.txt", sep=",", names=["energy", "trans"])


plt.figure()
plt.plot(trans_mg["energy"], trans_mg["trans"], label="Mg", color="blue")
plt.plot(trans_te["energy"], trans_te["trans"], label="Te", color="orange")
plt.xlabel("Energy (eV)")
plt.ylabel("Transmission")
plt.grid()
plt.legend()
plt.savefig('images/images_atto/calculated_Mg_Te_trans.png')
plt.show()


# =============================================================================
# 15
# =============================================================================

def model_f(x, w_in, x0, A, B):
    return A*np.exp(-2*(x-x0)**2/(w_in**2))+B 

norm_Y_no=sumY_no/max(sumY_no)



arg_sumY_no=np.arange(len(sumY_no))
popts, pcovs=curve_fit(model_f, arg_sumY_no, norm_Y_no, p0=[221, 350, 1, 0])
w_opts, x0_opts, A_opts, B_opts=popts


y_model_s=model_f(arg_sumY_no, w_opts, x0_opts, A_opts, B_opts)

calib=arg_sumY_no*40/1216   #in mm 
w2=w_opts*40/1215           #in mm 

plt.figure()
plt.plot(calib, norm_Y_no, label="sum along horizon. axis")
plt.plot(calib, y_model_s, label="fitting")
plt.xlabel("Calibrated distance, (mm)")
plt.ylabel("Intensity (arb. units)")
plt.legend()
plt.grid()
plt.savefig('images/images_atto/Spatial_inten_vs_calib_dis.png')
plt.show()
# =============================================================================
# 16
# =============================================================================
D=2500   #in mm 
Theta_target=w2/D

N1=w2/w_in

Theta_XUV=N1**2*Theta_target

# =============================================================================
# 17
# =============================================================================

Theta_IR=w_in/D
f=np.pi*w2*w_in/lmbda

# =============================================================================
# Calibration of pexel to energy eV 
# =============================================================================

E_start=42 # in eV 1012
d_E=2.4    # in eV 

peaks=[866, 915, 958, 1012, 1072, 1138, 1212, 1297]
h_start=35

energies=[44.4]
H_No=[37]
for i in range(7):
    e=E_start-i*d_E
    h=h_start-i*2
    energies.append(e)
    H_No.append(h)
    
# peaks_Y=[]
# Anno_data = {'H_order': H_No,
#              'peaks_Y': peaks_Y,
#         'x_coor': peaks}

# max_frame = pd.DataFrame(Anno_data)


# plt.figure(1)
# for i, row in max_frame.iterrows():
#     plt.annotate(f'{row["H_order"]:.2f}',
#                  xy=(row['x_coor'], row['peaks_Y']),
#                  xytext=(55, 1), textcoords='offset points', ha='right', va='bottom', color="green")
    

def model_curve(x, a, b):
    return a/(x)+b

popt_curve, pcov_curve=curve_fit(model_curve, peaks, energies, p0=[1, 0.05])
a_opt, b_opt=popt_curve

x_model_curve=np.linspace(min(peaks)-100, max(peaks)+200, 500)
y_model_curve=model_curve(x_model_curve, a_opt, b_opt)

plt.figure()
plt.scatter(peaks, energies, marker="x", color="red", label="peaks")
plt.plot(x_model_curve, y_model_curve, label="fitting")
plt.legend()
plt.grid()
plt.xlabel('pixel')
plt.ylabel('Energy, eV')
plt.legend()
plt.savefig('images/images_atto/Energy_pixel_calib.png')
plt.show()

def calib(x):
    return model_curve(x, a_opt, b_opt)


mg_clean_ener=trans_mg["energy"]
te_clean_ener=trans_te["energy"]

mg_clean_trans=trans_mg["trans"]/4
te_clean_trans=trans_te["trans"]


t1=np.argmin(mg_clean_ener<30)
t2=np.argmin(mg_clean_ener<64)

m1=np.argmin(te_clean_ener<30)
m2=np.argmin(te_clean_ener<64)

fig, ax1 = plt.subplots(figsize=(9, 6))
ax1.set_xlabel('Energy eV')
ax1.set_ylabel('Spectra', color="purple")
ax1.plot(calib(x[622:1312]), sumX_no[622:1312]/max( sumX_no[622:1312]), label="no filter*", color='purple', linewidth=1)
# ax1.plot(calib(x[622:1312]), sumX_mg[622:1312]/(2*max(sumX_mg[622:1312])), label="mg filter", color='black', linewidth=1)
ax1.plot(calib(x[622:1312]), sumX_te[622:1312]/max(sumX_te[622:1312]), label="te filter", linewidth=1)
ax1.tick_params(axis='y', labelcolor="purple")
ax1.legend(loc='upper left')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Transmission', color="blue")  # we already handled the x-label with ax1
# ax2.plot(calib(x_no), np.abs(T_w), label="Mg trans. Ex.", color="blue", linewidth=3)
ax2.plot(calib(x_no), np.abs(T_te), label="Te trans. Ex.", color="blue", linewidth=2)
# ax2.plot(mg_clean_ener[t1:t2], mg_clean_trans[t1:t2], linestyle='dashed', label='Mg Trans. cal', color="blue", linewidth=3)
ax2.plot(te_clean_ener[m1:m2], te_clean_trans[m1:m2], linestyle='dashed', label='Te Trans, cal', color="black", linewidth=3)
ax2.plot()
ax2.tick_params(axis='y', labelcolor="blue")
ax2.legend(loc='upper right')
plt.savefig('images/images_atto/all_together_final.png')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


# plt.figure()
# plt.subplot(211)
# plt.plot(calib(x[622:1312]), sumX_no[622:1312]/max( sumX_no[622:1312]), label="no filter*")
# plt.plot(calib(x[622:1312]), sumX_mg[622:1312]/max(sumX_mg[622:1312]), label="mg filter")
# plt.plot(calib(x[622:1312]), sumX_te[622:1312]/max(sumX_te[622:1312]), label="te filter")
# plt.xlabel("Energy (eV)")
# plt.ylabel("Intensity normalized ")
# plt.legend()
# plt.grid()  

# plt.subplot(212)
# plt.plot(calib(x_no), np.abs(T_w), label="Mg", color="blue")
# plt.plot(calib(x_no), np.abs(T_te), label="Te", color="orange")
# plt.xlabel("Energy (eV)")
# plt.ylabel("Transmission")
# plt.legend()
# plt.grid()
# plt.tight_layout()
# plt.savefig('images/images_atto/Calib_final.png')
# plt.show()



    
