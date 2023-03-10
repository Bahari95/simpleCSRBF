## coding
# just RBF using MultiQuadric
# Copyright 2023 BAHARI Mustapha

import matplotlib.pyplot            as     plt
from   mpl_toolkits.axes_grid1      import make_axes_locatable
from   mpl_toolkits.mplot3d         import axes3d
from   matplotlib                   import cm
from   mpl_toolkits.mplot3d.axes3d  import get_test_data
from   matplotlib.ticker            import LinearLocator, FormatStrFormatter
#..
from   scipy.sparse                 import kron
from   scipy.sparse                 import csr_matrix
from   scipy.sparse                 import csc_matrix, linalg as sla
from   numpy                        import zeros, linalg, asarray
from   numpy                        import cos, sin, pi, exp, sqrt, arctan2
from   tabulate                     import tabulate
from numpy                          import loadtxt 
import numpy                        as     np
import time

#./// TOOLS for CSRBF simulation
from gallery_section_02             import assemble_rbf_stiffnes_rhs
from simple_CSRBF                   import results

# le temps maximal
t_max     = 3.5 
x_max     = 2.
# The RBF function coefficient : scaling parameter
c         = 0.05

# ... discrètisation
coef_diff = 0.13/(0.11*7.5)
ro        = 0.15
#..
dx        = 0.1 
dt        = (ro*dx*dx)/coef_diff
nt        = int(t_max/dt)
nx        = int(x_max/dx)
t         = np.linspace(0, t_max, nt)
x         = np.linspace(0, x_max, nx)
T, X      = np.meshgrid(t, x)       


# stiffness matrix
stiffness  = np.zeros((nx,nt,nx,nt), dtype = np.double)   
# right hand side
rhs        = np.zeros((nx,nt), dtype = np.double)  

# ... Assembles matrix and rhs of RBF-Poisson
assemble_rbf_stiffnes_rhs(X, T, nx, nt, c, coef_diff, stiffness, rhs)

# ... Linear system from RBF 
stiffness = stiffness.reshape(nx*nt, nx*nt)

# ...
rhs       = rhs.reshape(nx*nt)

# ... Resolution of linear system
lu        = sla.splu(csc_matrix(stiffness))

# ...
alpha     = lu.solve(rhs)


# ... Computation of the RBF approximate solution
U_CSRBF_ET, u_exact = results(X, T, nx, nt, c, alpha, MQ = True, u_exact =  True)

         
print(" ERROR INFTY =", np.max(np.absolute(U_CSRBF_ET- u_exact)) )

##  plot
figtitle  = 'RBF TIME-SPACE FOR FOURIER EQUATION'

fig, axes = plt.subplots( 1, 2, figsize=[12,12], gridspec_kw={'width_ratios': [2.75, 2]} , num=figtitle )
for ax in axes:
   ax.set_aspect('equal')

#axes[0].set_title( 'Physical domain ' )
for i in range(nt):
    phidx = X[:,i]
    phidy = T[:,i]

    axes[0].plot(phidx, phidy, '-b', linewidth = 0.75)
for i in range(nx):
    phidx = X[i,:]
    phidy = T[i,:]

    axes[0].plot(phidx, phidy, '-b', linewidth = 0.75)
#axes[0].axis('off')
axes[0].margins(0,0)
im    = axes[1].contourf( X, T, U_CSRBF_ET, cmap= 'jet')
divider = make_axes_locatable(axes[1]) 
cax   = divider.append_axes("right", size="5%", pad=0.05, aspect = 40) 
plt.colorbar(im, cax=cax)
fig.tight_layout()
plt.subplots_adjust(wspace=0.3)
plt.savefig('r_refinement_ex.png')
plt.show()

# set up a figure twice as wide as it is tall
fig = plt.figure(figsize=plt.figaspect(0.5))
#===============
# First subplot
# set up the axes for the first plot
ax = fig.add_subplot(1, 2, 1, projection='3d')
# plot a 3D surface like in the example mplot3d/surface3d_demo
surf0 = ax.plot_surface(X, T, U_CSRBF_ET, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
#ax.set_xlim(0.0, 1.0)
#ax.set_ylim(0.0, 1.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#ax.set_title('Approximate solution in uniform mesh')
ax.set_xlabel('X',  fontweight ='bold')
ax.set_ylabel('Y',  fontweight ='bold')
# Add a color bar which maps values to colors.
fig.colorbar(surf0, shrink=0.5, aspect=25)

#===============
# Second subplot
ax = fig.add_subplot(1, 2, 2, projection='3d')
surf = ax.plot_surface(X, T, u_exact, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
#ax.set_xlim(0.0, 1.0)
#ax.set_ylim(0.0, 1.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#ax.set_title('Approximate Solution in adaptive meshes')
ax.set_xlabel('F1',  fontweight ='bold')
ax.set_ylabel('F2',  fontweight ='bold')
fig.colorbar(surf, shrink=0.5, aspect=25)
#plt.savefig('Poisson3D.png')
plt.show()
