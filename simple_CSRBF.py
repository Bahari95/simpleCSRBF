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

@types('real[:,:]', 'real[:,:]', 'int', 'int', 'real', 'int[:,:,:,:]', 'int[:,:]', 'real[:,:,:]')
def CSRBF_tools(X, Y, nx, ny, s, span, support, r_span):
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
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i-k1,j+k2])**2 + (Y[i,j]-Y[i-k1,j+k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span] = i-k1
	                 span[i,j,1,ij_span] = j+k2
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1

	    for k1 in range(0,nx-i):
	        for k2 in range(1,j+1):
	            r = sqrt((X[i,j]-X[i+k1,j-k2])**2 + (Y[i,j]-Y[i+k1,j-k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span]  = i+k1
	                 span[i,j,1,ij_span]  = j-k2
	                 r_span[i,j,ij_span]  = r
	                 ij_span             += 1
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i+k1,j+k2])**2 + (Y[i,j]-Y[i+k1,j+k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span] = i+k1
	                 span[i,j,1,ij_span] = j+k2
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1
	    support[i, j] = ij_span
	return 0.
# ...
pyccel_CSRBF_tools = epyccel(CSRBF_tools)

def CSRBF_basis(X, Y, nx, ny, s):
	# X is x-coordinates 2D matrix
	# Y is y-coordinates 2D matrix
	r_xy    = zeros((nx,ny, nx*ny))                           # ... computed distance for non niglicted RBF
	span    = zeros((nx,ny, 2, nx*ny), dtype = int)           # ... Determine the global index  for non-vanishing CSRBF       
	support = zeros((nx,ny), dtype = int)                     # ... Determine the span index for spans
	
	#...
	pyccel_CSRBF_tools(X, Y, nx, ny, s, span, support, r_xy)
	
	#---------------------------------
	#... resize variables
	max_span         = np.max(support)
	op_r_xy          = zeros((nx,ny, max_span))                 
	op_r_xy[:,:,:]   = r_xy[:,:,:max_span]
	del r_xy
	op_span          = np.zeros((nx,ny,2,max_span), dtype = int)     
	op_span[:,:,:]   = span[:,:,:,:max_span]
	del span
	return op_span, op_r_xy, support
	
	
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
     
assemble_mass_matrix = epyccel(assemble_mass)    

def results(X, Y, nx, ny, s, span, r_xy, support, control_points):
	#... Mass matrix 
	K          = np.zeros((nx,ny,nx,ny), dtype = np.double)

	# ... Assembles mass matrix for compute a soluotion in grid points
	assemble_mass_matrix(nx, ny, s, span, r_xy, support, K)

	K         = K.reshape(nx*ny, nx*ny)
	# ...
	u_csrbf   = K.dot(control_points)
	u_csrbf   = u_csrbf.reshape(nx,ny)
	# ...
	return u_csrbf



