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
from gallery_section_01             import assemble_Poisson_stiffnes_rhs_1D
from simplemeshless                 import CSRBF_basis_1D
from simplemeshless                 import results_1D
from simplemeshless                 import StencilMatrix_1D

#======================================================================================
##                    II- CSRBF method -II
## ... Spectrum in space: we use only the points (base) whose distance is less than s,
## in order to obtain the stencil matrix.
## ...see Wu and Wendland
#======================================================================================
# ..           Fourier-equation  : -\laplace(u) = f
#======================================================================================

# _ We use the Wendland function of the following form
# _f_ro(r)= (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2) où r=norm(X-Y,2)._
#.. s is a support
#--------------------------------------------------------------------------------------

#.. The number of points in x direction
nx    = 200             

# uniform cartesian mesh
X    = np.linspace(0,1, nx)

# ... control Support radius
# The RBF function coefficient : scaling parameter
s          = 0.15

# ... Computation of CSRBF TOOLS
max_span, span, r_x, support = CSRBF_basis_1D(X, nx, s)

# stiffness matrix
stiffness  = StencilMatrix_1D(nx, max_span, support, span)
# right hand side
rhs        = np.zeros(nx, dtype = np.double) 

# ... Assembles matrix and rhs of RBF-Poisson
assemble_Poisson_stiffnes_rhs_1D(X, nx, s, span, r_x, support, stiffness._data, rhs)

# ... Linear system from RBF 
stiffness = stiffness.tosparse()

# visualize the sparse matrix with Spy
plt.spy(stiffness)
plt.show()

# ... Resolution of linear system
lu        = sla.splu(csc_matrix(stiffness))
# ...
U         = lu.solve(rhs)

# ... Computation of the RBF approximate solution
u_csrbf_ts = results_1D(X, nx, s, U, max_span = max_span, span = span, r_x = r_x, support = support)

#... test 0
#u_exact = exp(-200.*(((X-.5)/0.4)**2+((Y-.5)/0.4)**2-0.6)**2 )

#... test 1
u_exact = sin(pi*X)
         
print("npoints = {} scaling parameter = {} max_span = {} Least square error ={}".format( nx, s, max_span, np.sqrt(np.sum(( u_csrbf_ts - u_exact)**2)) ) )

##  plot
plt.figure()
plt.axes().set_aspect('equal')
plt.plot( X, u_csrbf_ts,  color='r', lw = 2.5, ls='-', marker='s', markersize=1, markerfacecolor = "k", label='$\mathbf{CSRBF}$')
plt.plot( X, u_exact, color='b', lw = 2.5, ls='-', marker='v', markersize=1, markerfacecolor = "k", label = '$\mathbf{EXACT}$')
plt.xlabel('X',  fontweight ='bold')
plt.ylabel('f(X)',  fontweight ='bold')
plt.grid(color='b', linestyle='--', linewidth=0.5)
plt.margins(0.02,0.02)
plt.legend()
plt.show()
