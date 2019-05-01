
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties



r0r1 = (np.loadtxt('r.dat'))
r0   = r0r1[:,0]
r1   = r0r1[:,1]

xq1 = (np.loadtxt('q1.dat'))
x1   = xq1[:,0]
q1  = xq1[:,1]


xq2 = (np.loadtxt('q2.dat'))
x2   = xq2[:,0]
q2  = xq2[:,1]


xq3 = (np.loadtxt('q2.dat'))
x3   = xq3[:,0]
q3  = xq3[:,1]


plt.plot(r0,r1,label='Numeric r1')

plt.plot(x1,q1,label='Numeric q1')
plt.plot(x2,q2,label='Numeric q2')
plt.plot(x3,q3,label='Numeric q3')
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