# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 10:39:57 2023

@author: Owner
"""

import numpy as np 
import matplotlib.pyplot as plt 

tp=20   # fs 
lmbda=800   # nm 
c=299.73    # nm/fs

T=lmbda/c

N=2**18

dt=T/N 

t=np.arange(-T/2, T/2, dt) 
w=2*np.pi/T 

Et=np.cos(w*t)*np.exp(-(t/tp)**2+1j*w*t)

plt.figure()
plt.plot(t, Et)
plt.xlabel("time")
plt.ylabel("Electrical field")
plt.grid()
plt.savefig("Homework_9_selected_topics.png")
plt.show()