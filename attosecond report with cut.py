# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 18:19:22 2023

@author: Owner
"""

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit


plt.close('all')

image_path = 'Data/Applications_atto/Task11_holey_mirror_cut.png'
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

# def outliers(array):
#     mod_array=array.copy()
#     mod_array[1244:1488]=np.nan       
#     return mod_array

# def outliers2(array):
#     mod_array=array.copy()
#     mod_array[531:737]=np.nan       
#     return mod_array

# norm_x_clean=outliers(norm_x)
# norm_y_clean=outliers2(norm_y)

# # curve fitting: 
# # define the model function: 
# def model_f(x, w_in, x0, A, B):
#     return A*np.exp(-2*(x-x0)**2/(w_in**2))+B 

# valid_indices = ~np.isnan(norm_x_clean)

# # curve fit: 
# popt, pcov=curve_fit(model_f, len_x[valid_indices], norm_x[valid_indices], p0=[221, 1378, 1, 0])
# w_opt, x0_opt, A_opt, B_opt=popt

# x_model=np.arange(min(len_x), max(len_x), 1)
# y_model=model_f(x_model, w_opt, x0_opt, A_opt, B_opt)

# valid_indices2 = ~np.isnan(norm_y_clean)

# popt2, pcov2=curve_fit(model_f, len_y[valid_indices2], norm_y[valid_indices2], p0=[221, 500, 1, 0])
# w_opt2, x0_opt2, A_opt2, B_opt2=popt2


# x_model2=np.arange(min(len_y), max(len_y), 1)
# y_model2=model_f(x_model2, w_opt2, x0_opt2, A_opt2, B_opt2)

plt.figure()
plt.subplot(211)
plt.scatter(len_x, norm_x)

plt.subplot(212)
plt.scatter(len_y, norm_y)
plt.show()

# plt.figure()
# plt.subplot(211)
# plt.scatter(len_x, norm_x_clean)
# plt.plot(x_model, y_model, "b--")

# plt.subplot(212)
# plt.scatter(len_y, norm_y_clean)
# plt.plot(x_model2, y_model2, "r--")
# plt.show()

