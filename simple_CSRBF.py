## coding
#
# Copyright 2023 BAHARI Mustapha

"""
Basic module that provides the means for evaluating the CSRBF basis
functions. In order to simplify automatic Fortran code
generation with Pyccel.
"""
__all__ = ['CSRBF_tools']

from pyccel.decorators              import types
from pyccel.epyccel                 import epyccel
from numpy                          import zeros
import numpy as np

@types('real[:,:]', 'real[:,:]', 'int', 'int', 'real', 'int[:,:,:]', 'int[:,:,:]', 'int[:,:]', 'real[:,:,:]')
def CSRBF_tools(X, Y, nx, ny, s, span_x, span_y, spect_r, r_span):
	from numpy import sqrt
	for i in range(nx):
	  for j in range(ny):
	    ij_span = 0
	    for k1 in range(1,i+1):
	        for k2 in range(1,j+1):
	            r = sqrt((X[i,j]-X[i-k1,j-k2])**2 + (Y[i,j]-Y[i-k1,j-k2])**2)
	            if r < s :
	                 span_x[i,j,ij_span] = i-k1
	                 span_y[i,j,ij_span] = j-k2
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i-k1,j+k2])**2 + (Y[i,j]-Y[i-k1,j+k2])**2)
	            if r < s :
	                 span_x[i,j,ij_span] = i-k1
	                 span_y[i,j,ij_span] = j+k2
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1

	    for k1 in range(0,nx-i):
	        for k2 in range(1,j+1):
	            r = sqrt((X[i,j]-X[i+k1,j-k2])**2 + (Y[i,j]-Y[i+k1,j-k2])**2)
	            if r < s :
	                 span_x[i,j,ij_span]  = i+k1
	                 span_y[i,j,ij_span]  = j-k2
	                 r_span[i,j,ij_span]  = r
	                 ij_span             += 1
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i+k1,j+k2])**2 + (Y[i,j]-Y[i+k1,j+k2])**2)
	            if r < s :
	                 span_x[i,j,ij_span] = i+k1
	                 span_y[i,j,ij_span] = j+k2
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1
	    spect_r[i, j] = ij_span
	return 0.
# ...
pyccel_CSRBF_tools = epyccel(CSRBF_tools)

def CSRBF_basis(X, Y, nx, ny, s):
	# X is x-coordinates 2D matrix
	# Y is y-coordinates 2D matrix
	r_xy    = zeros((nx,ny, nx*ny))                        # ... computed distance for non niglicted RBF
	span_x  = zeros((nx,ny, nx*ny), dtype = int)           # ... Determine the global index for x coordinate for non-vanishing CSRBF       
	span_y  = zeros((nx,ny, nx*ny), dtype = int)           # ... Determine the global index for y coordinate for non-vanishing CSRBF
	spect_r = zeros((nx,ny), dtype = int)                  # ... Determine the span index for spans
	
	#...
	pyccel_CSRBF_tools(X, Y, nx, ny, s, span_x, span_y, spect_r, r_xy)
	
	#---------------------------------
	#... resize variables
	max_span         = np.max(spect_r)
	op_r_xy          = zeros((nx,ny, max_span))                 
	op_r_xy[:,:,:]   = r_xy[:,:,:max_span]
	del r_xy
	op_span_x        = np.zeros((nx,ny, max_span), dtype = int)     
	op_span_x[:,:,:] = span_x[:,:,:max_span]
	del span_x
	op_span_y        = np.zeros((nx,ny, max_span), dtype = int)      
	op_span_y[:,:,:] = span_y[:,:,:max_span]
	del span_y
	return op_span_x, op_span_y, op_r_xy, spect_r
