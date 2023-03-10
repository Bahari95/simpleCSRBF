## coding
# CSRBF USING Wendland function C6
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
from gallery_section_02             import assemble_csrbf_stiffnes_rhs
from simplemeshless                 import CSRBF_basis
from simplemeshless                 import results
from simplemeshless                 import StencilMatrix


#======================================================================================
##                    II- CSRBF method -II
## ... Spectrum in space: we use only the points (base) whose distance is less than s,
## in order to obtain the stencil matrix.
## ...see Wu and Wendland
#======================================================================================
# ..           Fourier-equation  : \partial_t u - coefs \partial_x(u)=f
#======================================================================================

# le temps maximal
t_max    = 2.5    
x_max    = 2.
# The RBF function coefficient : scaling parameter
s         = 0.3

# ... discrètisation
coef_diff = 0.1576
ro        = 0.15
#..
#dx        = 0.05 
#dt        = 0.05 #(ro*dx*dx)/coef_diff
nt        = 100 #int(t_max/dt)
nx        = 100 #int(x_max/dx)
t         = np.linspace(0, t_max, nt)
x         = np.linspace(0, x_max, nx)
T, X      = np.meshgrid(t, x)

# ... Computation of CSRBF TOOLS
max_span_x, max_span_t, span, r_xt, support = CSRBF_basis(X, T, nx, nt, s)

# stiffness matrix
stiffness  = StencilMatrix(nx, nt, max_span_x, max_span_t, support, span)

# right hand side
rhs        = np.zeros((nx,nt), dtype = np.double) 

# ... Assembles matrix and rhs of RBF-Poisson
assemble_csrbf_stiffnes_rhs(X, T, nx, nt, s, span, r_xt, support, coef_diff, stiffness._data, rhs)

# ... Linear system from RBF 
stiffness = stiffness.tosparse()
#stiffness = stiffness.reshape(nx*ny, nx*ny)

# visualize the sparse matrix with Spy
#plt.spy(stiffness)
#plt.show()
# ...
rhs       = rhs.reshape(nx*nt)

# ... Resolution of linear system
lu        = sla.splu(csc_matrix(stiffness))

# ...
U         = lu.solve(rhs)

# ... Computation of the RBF approximate solution
U_CSRBF_ET, u_exact = results(X, T, nx, nt, s, U, max_span_x = max_span_x, max_span_y = max_span_t, span = span, r_xy = r_xt, support = support, u_exact =  True)

         
print("npoints = {} scaling parameter = {} max_span = {} Least square error ={}".format(nx*nt, s, max(max_span_x,max_span_t), np.sqrt(np.sum(( U_CSRBF_ET - u_exact)**2)) ) )

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
#plt.savefig('Fourier-equation.png')
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
