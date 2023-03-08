## coding
#
# Copyright 2023 BAHARI Mustapha

"""
Basic module that provides the means for evaluating the CSRBF basis
functions. In order to simplify automatic Fortran code
generation with Pyccel.
"""
__all__ = ['CSRBF_tools',
          'assemble_mass']

from pyccel.decorators              import types
from pyccel.epyccel                 import epyccel
from numpy                          import zeros
import numpy as np

# -------
@types('real[:,:]', 'real[:,:]', 'int', 'int', 'real', 'int[:,:]')
def fin_support(X, Y, nx, ny, s, support):
	from numpy import sqrt
	for i in range(nx):
	  for j in range(ny):
	    ij_span = 0
	    for k1 in range(1,i+1):
	        for k2 in range(1,j+1):
	            r = sqrt((X[i,j]-X[i-k1,j-k2])**2 + (Y[i,j]-Y[i-k1,j-k2])**2)
	            if r < s :
	                 ij_span            += 1
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i-k1,j+k2])**2 + (Y[i,j]-Y[i-k1,j+k2])**2)
	            if r < s :
	                 ij_span            += 1

	    for k1 in range(0,nx-i):
	        for k2 in range(1,j+1):
	            r = sqrt((X[i,j]-X[i+k1,j-k2])**2 + (Y[i,j]-Y[i+k1,j-k2])**2)
	            if r < s :
	                 ij_span             += 1
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i+k1,j+k2])**2 + (Y[i,j]-Y[i+k1,j+k2])**2)
	            if r < s :
	                 ij_span            += 1
	    support[i, j] = ij_span
	return 0.
pyccel_fin_support = epyccel(fin_support)

# -------
@types('real[:,:]', 'real[:,:]', 'int', 'int', 'real', 'int', 'int[:,:,:,:]', 'real[:,:,:]')
def CSRBF_tools(X, Y, nx, ny, s, max_span, span, r_span):
	from numpy import sqrt
	for i in range(nx):
	  for j in range(ny):
	    
	    ij_span = 0
	    for k1 in range(1,i+1):
	        for k2 in range(1,j+1):
	            r = sqrt((X[i,j]-X[i-k1,j-k2])**2 + (Y[i,j]-Y[i-k1,j-k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span] = i-k1
	                 span[i,j,1,ij_span] = j-k2
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1
	            if ij_span == max_span:
	                    break
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i-k1,j+k2])**2 + (Y[i,j]-Y[i-k1,j+k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span] = i-k1
	                 span[i,j,1,ij_span] = j+k2
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1
	            if ij_span == max_span:
	                    break
	        if ij_span == max_span:
	               break
	    for k1 in range(0,nx-i):
	        for k2 in range(1,j+1):
	            r = sqrt((X[i,j]-X[i+k1,j-k2])**2 + (Y[i,j]-Y[i+k1,j-k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span]  = i+k1
	                 span[i,j,1,ij_span]  = j-k2
	                 r_span[i,j,ij_span]  = r
	                 ij_span            += 1
	            if ij_span == max_span:
	                    break
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i+k1,j+k2])**2 + (Y[i,j]-Y[i+k1,j+k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span] = i+k1
	                 span[i,j,1,ij_span] = j+k2
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1
	            if ij_span == max_span:
	                    break
	        if ij_span == max_span:
	               break
	return 0.
# ...
pyccel_CSRBF_tools = epyccel(CSRBF_tools)

def CSRBF_basis(X, Y, nx, ny, s):
	# X is x-coordinates 2D matrix
	# Y is y-coordinates 2D matrix
	support       = zeros((nx,ny), dtype = int)                  # ... Determine the span index for spans
	pyccel_fin_support(X, Y, nx, ny, s, support)
	# ...
	max_span      = int(np.max(support))
	#...
	r_xy          = zeros((nx,ny, max_span))                     # ... computed distance for non niglicted RBF                 
	span          = np.zeros((nx,ny,2,max_span), dtype = int)    # ... Determine the global index  for non-vanishing CSRBF
	#...
	pyccel_CSRBF_tools(X, Y, nx, ny, s, max_span, span, r_xy)

	return span, r_xy, support
	
# .... Wendland function C6
@types('int', 'int', 'real', 'int[:,:,:,:]', 'real[:,:,:]', 'int[:,:]', 'real[:,:,:,:]',)
def assemble_mass(nx, ny, s, span, r_xy, support, K): 
     from numpy import exp
     from numpy import cos
     from numpy import sin
     from numpy import pi
     from numpy import sqrt
     for i1 in range(0,nx):
        for i2 in range(0,ny):
                 
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1 = span[i1, i2, 0, ij_span]
                         j2 = span[i1, i2, 1, ij_span]
                         r  = r_xy  [i1, i2, ij_span]
                         #...
                         K[i1,i2,j1,j2]  = (1-r/s)**6.*(3+18.*r/s+35.*(r/s)**2)
     return 0
     
assemble_mass_matrixWC6 = epyccel(assemble_mass)    

# ..... MultiQuadric	
@types('int', 'int', 'real', 'int[:,:,:,:]', 'real[:,:,:]', 'int[:,:]', 'real[:,:,:,:]')
def assemble_massMQ(nx, ny, c, span, r_xy, support, K): 
     from numpy import exp
     from numpy import cos
     from numpy import sin
     from numpy import pi
     from numpy import sqrt
     for i1 in range(0,nx):
        for i2 in range(0,ny):
                 
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1 = span[i1, i2, 0, ij_span]
                         j2 = span[i1, i2, 1, ij_span]
                         r  = r_xy  [i1, i2, ij_span]
                         #...
                         K[i1,i2,j1,j2]  = sqrt(r**2+c**2)
     return 0
     
assemble_mass_matrixMQ = epyccel(assemble_massMQ)    

@types('real[:,:]', 'real[:,:]', 'int', 'int', 'real[:,:]')
def assemble_sol_exact(X_cor, Y_cor, nx, ny, u_exact): 
     from numpy import exp
     from numpy import cos
     from numpy import sin
     from numpy import pi
     from numpy import sqrt
     for i1 in range(0,nx):
        for i2 in range(0,ny):
                 
                 x  = X_cor[i1,i2]
                 t  = Y_cor[i1,i2] 
                 for n in range(1000):
                         u_exact[i1,i2] += 800./(pi**2*(2.*n+1)**2)*cos(pi*(2.*n-1)*(x-1.)*0.5)*exp(-0.3738*(2.*n+1)**2*t)
     return 0
sol_exact = epyccel(assemble_sol_exact) 

def results(X, Y, nx, ny, s, span, r_xy, support, control_points, MQ = None, u_exact = None):
	#... Mass matrix 
	K          = np.zeros((nx,ny,nx,ny), dtype = np.double)
	
	# ... Assembles mass matrix for compute a soluotion in grid points
	if MQ == None:
	       assemble_mass_matrixWC6(nx, ny, s, span, r_xy, support, K)
	else :
	       assemble_mass_matrixMQ(nx, ny, s, span, r_xy, support, K)
	K         = K.reshape(nx*ny, nx*ny)
	# ...
	u_csrbf   = K.dot(control_points)
	u_csrbf   = u_csrbf.reshape(nx,ny)
	# ...
	if u_exact is None :
	      return u_csrbf
	else :
	      u_exact  = np.zeros((nx,ny))
	      sol_exact(X, Y, nx, ny, u_exact)
	      return u_csrbf, u_exact 
