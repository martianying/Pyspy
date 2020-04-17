import numpy as np
import matplotlib.pyplot as plt
import math as m
import random 
from matplotlib.pyplot import figure

maxp = 5 * 10**3; 
scale = 3;
rcyl = 15;
rcylin = 0;
zmax = 4;
zmax5 = zmax * 0.5;

r_basket = []
x_basket = []
y_basket = []

i= 0

while i < maxp:
    
    riseed = random.random()
    ri = - scale * m.log( ( 1. - riseed) + m.exp( - rcyl/scale) * riseed) 
    
    if rcyl > ri and ri > rcylin:
        i = i + 1
        itheta = 2. * m.pi * random.random()
        xi = ri * m.cos(itheta)
        x_basket = x_basket + [xi]
        
        yi = ri * m.sin(itheta)
        y_basket = y_basket + [yi]
        
#         zi = 2.0 * ( zmax * ran2(iseed) - zmax5)
        r_basket = r_basket + [ri]
        
    else:
        continue
        
figure(num=None, figsize=(4, 4), dpi=300, facecolor='w', edgecolor='k')
plt.xlim(-20, 20)
plt.ylim(-20, 20)
plt.scatter(x_basket, y_basket, s = 0.2, color = 'black')    
        
    
    