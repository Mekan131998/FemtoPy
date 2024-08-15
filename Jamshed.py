# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 20:00:44 2023

@author: Owner
"""


import numpy as np
import shutil 
import matplotlib.pyplot as plt




textfiles = ["p_143437a.txt",
"p_143712a.txt",
"p_143942a.txt",
"p_144212a.txt",
"p_144442a.txt",
"p_144712a.txt",
"p_144942a.txt",
"p_145212a.txt",
"p_145442a.txt",
"p_145712a.txt",
"p_145942a.txt",
"p_150212a.txt",
"p_150442a.txt",
"p_150712a.txt",
"p_150943a.txt",
"p_151213a.txt",
"p_151443a.txt",
"p_151713a.txt",
"p_151943a.txt",
"p_152213a.txt",
"p_152443a.txt",
"p_152713a.txt",
"p_152943a.txt",
"p_153213a.txt",
"p_153443a.txt",
"p_153713a.txt",
"p_153943a.txt",
"p_154213a.txt",
"p_154443a.txt",
"p_154713a.txt",
"p_154943a.txt",
"p_155213a.txt",
"p_155443a.txt",
"p_155713a.txt",
"p_155943a.txt",
"p_160213a.txt",
"p_160444a.txt"]

file = "for 532nm.txt"

with open(file, "w") as wfd: 
    for f in textfiles: 
        with open(f, "r") as fd: 
            shutil.copyfileobj(fd, wfd) 
            

col0 = np.loadtxt(file, delimiter=',' , usecols = 0 )
col1 = np.loadtxt(file, delimiter=',' , usecols = 1 )
col2 = np.loadtxt(file, delimiter=',' , usecols = 2 )
# col3 = np.loadtxt(file, delimiter=',' , usecols = 3 )
# col4 = np.loadtxt(file, delimiter=',' , usecols = 4 )
    
t=range(len(col0))
# plt.figure()
# # plt.plot(t,col1,'r', label='Wavelength 1064nm' )
# plt.plot(t,col2, 'g', label='Wavelength 532nm')
# # plt.plot(t,col3, 'b', label='Wavelength 355nm')
# # plt.plot(t,col4, 'y', label='Wavelength 266nm')
# plt.legend(loc='upper right')
# plt.xlabel('Input laser wavelength [nm] ')
# plt.ylabel('Intensity [a.u.] ')
# np.savetxt('Caliberation of 532nm.txt',t, delimiter='')
# print("Data saved to 'Signal for 532nm.txt'.")
# plt.show()


# def smooth_points(data, angle_threshold):
#     angles = np.arctan(np.diff(data))  # Calculate angles between consecutive points
#     idx = []
#     jump_points = []

#     for i in range(1, len(angles)):
#         angle_change = np.abs(angles[i] - angles[i - 1])
#         if angle_change > np.abs(angle_threshold):
#             jump_points.append(data[i])
#             idx.append(i)

#     return jump_points, idx

# # Example usage:

# threshold = 1.39
# result, idx= smooth_points(col2, threshold)
# t=np.array(t)
# result=np.array(result)
# idx=np.array(idx)

# print("Jump points:", result)
# print(t[idx])

from scipy.signal import argrelextrema


def find_local_extrema(data):
    maxima = argrelextrema(data, np.greater)[0]
    minima = argrelextrema(data, np.less)[0]
    return maxima, minima

maxima, minima = find_local_extrema(col2)

extrem_data=[col2[i] for i in maxima]

# Plotting the data and the local maxima and minima
plt.plot(col2, label='Data')
plt.scatter(maxima, extrem_data , color='red', label='Local Maxima')
plt.legend()
plt.show()

# plt.figure(2)
# plt.plot(t,col2, 'g', label='Wavelength 532nm')
# plt.scatter(t[idx], result)
# plt.legend(loc='upper right')
# plt.xlabel('Input laser wavelength [nm] ')
# plt.ylabel('Intensity [a.u.] ')
# np.savetxt('Caliberation of 532nm.txt',t, delimiter='')
# print("Data saved to 'Signal for 532nm.txt'.")
# plt.show()

def find_jump_points(data, threshold):
    jump_points = []
    idx1=[]
    
    for i in range(1, len(data)):
        difference = data[i] - data[i-1]
        if difference > threshold:
            jump_points.append(data[i])
            idx1.append(i)
    
    return jump_points, idx1

threshold1=1000

points, idx1=find_jump_points(extrem_data, threshold1)
points, idx1=np.array(points), np.array(idx1)


plt.figure(3)
# plt.plot(t,col2, 'g', label='Wavelength 532nm')
plt.scatter(t[idx1], points)
plt.legend(loc='upper right')
plt.xlabel('Input laser wavelength [nm] ')
plt.ylabel('Intensity [a.u.] ')
np.savetxt('Caliberation of 532nm.txt',t, delimiter='')
print("Data saved to 'Signal for 532nm.txt'.")
plt.show()



