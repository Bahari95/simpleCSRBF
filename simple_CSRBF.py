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

from request_csrbf                  import pyccel_fin_support
from request_csrbf                  import pyccel_CSRBF_tools
from request_csrbf                  import assemble_mass_matrixWC6
# ...
from gallery                        import assemble_mass_matrixMQ
from gallery                        import sol_exact
# ...
from numpy                          import zeros
import numpy as np

def CSRBF_basis(X, Y, nx, ny, s):
	'''
	CSRBF_basis : computes the span according to scaling parameter "s" and euclidian distance for all points in the span,
	
	where,
	
	# X  is x-coordinates 2D matrix
	# Y  is y-coordinates 2D matrix
	# nx is number of points in x direction
	# ny is number of points in y direction
	# s  is scaling parameter 
	TODO : LOC_SPAN
	-> max_span_x : 
	-> max_span_y :
	'''
	support       = zeros((nx,ny), dtype = int)                  # ... Determine the span index for spans
	pyccel_fin_support(X, Y, nx, ny, s, support)
	# ...
	max_span      = int(np.max(support))
	#...
	r_xy          = zeros((nx,ny, max_span))                     # ... computed distance for non niglicted RBF                 
	span          = np.zeros((nx,ny,4,max_span), dtype = int)    # ... Determine the global index  for non-vanishing CSRBF
	#...
	pyccel_CSRBF_tools(X, Y, nx, ny, s, support, span, r_xy)
	# ...
	max_span_x, max_span_y = np.max(span[:,:,2,:])+1, np.max(span[:,:,3,:])+1
	return max_span_x, max_span_y, span, r_xy, support
	

def results(X, Y, nx, ny, s, control_points, span = None, r_xy = None, support = None, MQ = None, u_exact = None):
	'''
	this function compute the approximate solution by multiplying mass matrix and control points
	# X  is x-coordinates 2D matrix
	# Y  is y-coordinates 2D matrix
	# nx is number of points in x direction
	# ny is number of points in y direction
	# s  is scaling parameter 
	# control_points is a control points associated to the approximate solution
	# span contains the index of local and global points contained in the circle induced by the radius "s"
	# r_xy contains the euclidian distance between closed points according ton span
	# support contains a total index of points in circle induced by "s"
	# MQ for multiquadratic approach
	# MQ is none for Wendland function C6
	'''
	#... Mass matrix 
	K          = np.zeros((nx,ny,nx,ny), dtype = np.double)
	
	# ... Assembles mass matrix for compute a soluotion in grid points
	if MQ is None:
	       assemble_mass_matrixWC6(nx, ny, s, span, r_xy, support, K)
	else :
	       assemble_mass_matrixMQ(nx, ny, s, X, Y, K)
	K         = K.reshape(nx*ny, nx*ny)
	# ...
	u   = K.dot(control_points)
	u   = u.reshape(nx,ny)
	# ...
	if u_exact is None :
	      return u
	else :
	      u_exact  = np.zeros((nx,ny))
	      sol_exact(X, Y, nx, ny, u_exact)
	      return u, u_exact 
