## coding
#
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
from   numpy                        import loadtxt 
import numpy                        as     np
import time

#./// TOOLS for CSRBF simulation
from gallery_section_01             import assemble_Poisson_stiffnes_rhs
from simple_CSRBF                   import CSRBF_basis
from simple_CSRBF                   import results
from linalg                         import StencilMatrix

#======================================================================================
##                    II- CSRBF method -II
## ... Spectrum in space: we use only the points (base) whose distance is less than s,
## in order to obtain the stencil matrix.
## ...see Wu and Wendland
#======================================================================================

# _ We use the Wendland function of the following form
# _f_ro(r)= (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2) où r=norm(X-Y,2)._
#.. s is a support
#--------------------------------------------------------------------------------------

#.. The number of points in x direction
nx    = 100             
#.. The number of points in y direction
ny    = 100
# ... The total number of points (Number of nodes)
n     = nx * ny

# uniform cartesian mesh
xl    = np.linspace(0,1, nx)
yl    = np.linspace(0,1, ny)
Y, X  = np.meshgrid(yl, xl)

#... Adapted mesh using nx = 100 times ny = 100 and test 0 only
#nx = 100; ny =100; n =nx*ny;
#X = loadtxt('X_ad_'+str(nx)+'.txt')
#Y = loadtxt('Y_ad_'+str(ny)+'.txt')

# ... control Support radius
# The RBF function coefficient : scaling parameter
s          = 0.1

# ... Computation of CSRBF TOOLS
max_span_x, max_span_y, span, r_xy, support = CSRBF_basis(X, Y, nx, ny, s)

# stiffness matrix
stiffness  = StencilMatrix(nx, ny, max_span_x, max_span_y, support, span)
# right hand side
rhs        = np.zeros((nx,ny), dtype = np.double) 

# ... Assembles matrix and rhs of RBF-Poisson
assemble_Poisson_stiffnes_rhs(X, Y, nx, ny, s, span, r_xy, support, stiffness._data, rhs)

# ... Linear system from RBF 
stiffness = stiffness.tosparse()

# visualize the sparse matrix with Spy
#plt.spy(stiffness)
#plt.show()

# ...
rhs       = rhs.reshape(nx*ny)

# ... Resolution of linear system
lu        = sla.splu(csc_matrix(stiffness))
# ...
U         = lu.solve(rhs)

# ... Computation of the RBF approximate solution
U_CSRBF_ET = results(X, Y, nx, ny, s, U, max_span_x = max_span_x, max_span_y = max_span_y, span = span, r_xy = r_xy, support = support)

#... test 0
#U_exact_ET = exp(-200.*(((X-.5)/0.4)**2+((Y-.5)/0.4)**2-0.6)**2 )

#... test 1
U_exact_ET = sin(pi*X)*sin(pi*Y)
         
print(" Least square error =", np.sqrt(np.sum(( U_CSRBF_ET - U_exact_ET)**2)) )

##  plot
figtitle  = 'RBF TIME-SPACE FOR FOURIER EQUATION'

fig, axes = plt.subplots( 1, 2, figsize=[12,12], gridspec_kw={'width_ratios': [2.75, 2]} , num=figtitle )
for ax in axes:
   ax.set_aspect('equal')

#axes[0].set_title( 'Physical domain ' )
for i in range(ny):
    phidx = X[:,i]
    phidy = Y[:,i]

    axes[0].plot(phidx, phidy, '-b', linewidth = 0.75)
for i in range(nx):
    phidx = X[i,:]
    phidy = Y[i,:]

    axes[0].plot(phidx, phidy, '-b', linewidth = 0.75)
#axes[0].axis('off')
axes[0].margins(0,0)
im    = axes[1].contourf( X, Y, U_CSRBF_ET, cmap= 'jet')
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
surf0 = ax.plot_surface(X, Y, U_exact_ET, rstride=1, cstride=1, cmap=cm.coolwarm,
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
surf = ax.plot_surface(X, Y, U_CSRBF_ET, cmap=cm.coolwarm,
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
