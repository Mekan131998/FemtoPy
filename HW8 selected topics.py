# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 17:20:34 2023

@author: Owner
"""

import numpy as np 
import matplotlib.pyplot as plt 

plt.close('all')
L=10
x=np.linspace(0, L, 1000)

def density(n):
    den=0
    for i in range(n+1):
        den+=(2/L)*(np.sin(i*np.pi*x/L))**2
    return den

def norm(f):
    return f/max(f)

n1=norm(density(1))
n2=norm(density(2))
n3=norm(density(10))
n4=norm(density(100))

plt.figure()
plt.subplot(221)
plt.plot(x, n1, "b", label="n=1")
plt.ylabel("Density (a.u)")
plt.xlabel("Distance x (a.u)")
plt.legend()
plt.grid()

plt.subplot(222)
plt.plot(x, n2, "b", label="n=2")
plt.ylabel("Density (a.u)")
plt.xlabel("Distance x (a.u)")
plt.legend()
plt.grid()

plt.subplot(223)
plt.plot(x, n3, "b", label="n=10")
plt.ylabel("Density (a.u)")
plt.xlabel("Distance x (a.u)")
plt.legend()
plt.grid()

plt.subplot(224)
plt.plot(x, n4, "b", label="n=100")
plt.ylabel("Density (a.u)")
plt.xlabel("Distance x (a.u)")
plt.legend()
plt.grid()
plt.savefig('images/HW8')
plt.tight_layout()