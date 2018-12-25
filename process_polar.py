import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from matplotlib import cm
import pandas as pd

u = pd.read_csv('u_mesh.dat',header=None, sep='\s+').values[:,0]
u_exact = pd.read_csv('u_mesh_exact.dat', header=None, sep='\s+').values[:,0]

N = len(u)

r = np.linspace(0,1,N+1) # generate with the boundary
r = r[0:-1]

plt.plot(r,u,'r--',label='expr')
plt.plot(r,u_exact,'k-',label='exact')
plt.legend(loc=1)

dnorm = np.sum( 2*np.pi*r/N*(u-u_exact)**2 )

plt.text(0.2,0.4,s=r'$|u-u_{exact}|_{L_2}$'+' = %e'%dnorm)

savename = str(input('give a save name: ')) + '.png'
plt.savefig(savename)
plt.show()
