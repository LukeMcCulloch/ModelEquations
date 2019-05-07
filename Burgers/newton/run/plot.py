
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


sclices = {1:[0,62],
           2:[62,124],
           3:[124,186],
           4:[186,248],
           5:[248,310],
           6:[310,372]}

xq1 = (np.loadtxt('EulerImplicit_FV.dat'))
x1   = xq1[:,0]
q1  = xq1[:,1]


for sl in sclices:
    rg = sclices[sl]
    plt.plot(x1[rg[0]:rg[1]],
             q1[rg[0]:rg[1]],label='Numeric U')
    plt.axis(xmin=-5., ymin= -5., xmax=20.0, ymax=5.)

plt.grid()


L2norm = (np.loadtxt('L2norm.dat'))
x   = L2norm[:,0]
L2  = L2norm[:,1]
plt.plot(x,L2,label='L2 norm')
plt.xlabel(" x position, m")
plt.ylabel("L2 norm")
plt.title('$L_{2}$ norm newton convergence')
font = FontProperties(stretch='condensed', size='small')
plt.grid()