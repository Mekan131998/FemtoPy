# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:10:08 2023

@author: Owner
"""
# SF10 
# 

# BK7
# 
# Fused silica 

import numpy as np 


c=299.792 

SF10=[1, 0.38, 2.5]
BK7=[2, 0.3, 2.5]
FS=[3, 0.21, 6.7]


def refrindex(w, glassType):
    if w==0: 
        return 1 
    
       elif glassType[0]==1:
            w_min=2*np.pi*c/glassType[-1]
            w_max=2*np.pi*c/glassType[-2]
            if (w>0 and w<w_min) and w>w_max:
                return 1
            else: 
                x=2*np.pi*c/(w*1000)
                n=(1+1.62153902/(1-0.0122241457/x**2)+0.256287842/(1-0.0595736775/x**2)+1.64447552/(1-147.468793/x**2))**.5
                return n            
       elif glassType[0]==2:
            w_min=2*np.pi*c/glassType[-1]
            w_max=2*np.pi*c/glassType[-2]
            if (w>0 and w<w_min) and w>w_max:
                return 1
            else: 
                x=2*np.pi*c/(w*1000)
                n=(1+1.03961212/(1-0.00600069867/x**2)+0.231792344/(1-0.0200179144/x**2)+1.01046945/(1-103.560653/x**2))**.5
                return n
        elif glassType[0]==3:
            w_min=2*np.pi*c/glassType[-1]
            w_max=2*np.pi*c/glassType[-2]
            if (w>0 and w<w_min) and w>w_max:
                return 1
            else: 
                x=2*np.pi*c/(w*1000)
                n=(1+0.6961663/(1-(0.0684043/x)**2)+0.4079426/(1-(0.1162414/x)**2)+0.8974794/(1-(9.896161/x)**2))**.5
            return n

lmbda=0.8
w=2*np.pi*c/lmbda
w_min=2*np.pi*c/FS[-1]
w_max=2*np.pi*c/FS[-2]
print("w_min", w_min)
print("w_max", w_max)
print("w=",w )
print("n=", refrindex(w, FS))

# increase the number of the array elements 
# decrease the width of the angular frequency window 

# TOD not so large GDD is dominant term 
# subtitude GDD value this value to the pulse duration formula 
# t_min=GD0-2*tau 
# t_max=GD0+2*tau