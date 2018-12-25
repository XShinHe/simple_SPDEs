import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from matplotlib import cm
import pandas as pd

u = pd.read_csv('u_mesh.dat',header=None, sep='\s+').values
u_exact = pd.read_csv('u_mesh_exact.dat', header=None, sep='\s+').values

N = len(u)

x,y = np.mgrid[0:1:(N*1j+1j),0:1:(N*1j+1j)] # generate with the boundary

x = x[0:-1,0:-1]
y = y[0:-1,0:-1]

ax = plt.subplot(111, projection='3d')
ax.set_title('Solution on first quadrant (cylindrically symmetric)');
#ax.plot_surface(x,p,rho,rstride=2, cstride=1, cmap=cm.jet)
ax.plot_surface(x,y,u,rstride=2, cstride=1, cmap=cm.jet)
ax.plot_surface(x,y,u_exact,rstride=2, cstride=1, cmap='gnuplot')
#ax.plot_surface(x,p,-rhoexact,rstride=2, cstride=1, cmap=cm.jet)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u')

savename = str(input('give a save name: ')) + '.png'
plt.savefig(savename)
plt.show()

