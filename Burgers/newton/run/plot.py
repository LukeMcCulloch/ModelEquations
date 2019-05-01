
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties




xq1 = (np.loadtxt('EulerImplicit_FV.dat'))
x1   = xq1[:,0]
q1  = xq1[:,1]



plt.plot(x1,q1,label='Numeric U')
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