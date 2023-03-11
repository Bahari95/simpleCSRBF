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

from request_csrbf                  import pyccel_fin_support_1D
from request_csrbf                  import pyccel_CSRBF_tools_1D
from request_csrbf                  import assemble_mass_matrixWC6_1D
from linalg                         import StencilMatrix_1D
# ..
from request_csrbf                  import pyccel_fin_support
from request_csrbf                  import pyccel_CSRBF_tools
from request_csrbf                  import assemble_mass_matrixWC6
from linalg                         import StencilMatrix
# ...
from gallery                        import assemble_mass_matrixMQ
from gallery                        import sol_exact

# ...
from numpy                          import zeros
import numpy as np

# --------------------------
# In one dimension case
# --------------------------
def CSRBF_basis_1D(X, nx, s):
	'''
	CSRBF_basis : computes the span according to scaling parameter "s" and euclidian distance for all points in the span,
	
	where,
	
	# X  is x-coordinates 1D matrix
	# nx is number of points
	# s  is scaling parameter 
	-> max_span : 
	'''
	support       = zeros(nx, dtype = int)               # ... Determine the span index for spans
	pyccel_fin_support_1D(X, nx, s, support)
	# ...
	max_span      = int(np.max(support))
	# ...
	r_x           = zeros((nx, max_span))                # ... computed distance for non niglicted RBF                 
	span          = zeros((nx, max_span), dtype = int)   # ... Determine the global index  for non-vanishing CSRBF
	# ...
	pyccel_CSRBF_tools_1D(X, nx, s, support, span, r_x)
	# ...
	return max_span, span, r_x, support

def results_1D(X, nx, s, control_points,  max_span = None, span = None, r_x = None, support = None):
	'''
	this function compute the approximate solution by multiplying mass matrix and control points
	# X  is x-coordinates 1D matrix
	# nx is number of points
	# s  is scaling parameter 
	# control_points is a control points associated to the approximate solution
	# max_span is the maximum number of points close to the objective point
	# span contains the index of local and global points contained in the circle induced by the radius "s"
	# r_x contains the euclidian distance between closed points according ton span
	# support contains a total index of points in circle induced by "s"
	# WC6 for Wendland function C6
	'''
	#... Mass matrix  sparse matrix TODO
	K   = StencilMatrix_1D(nx, max_span, support, span)
	
	# ... Assembles mass matrix for compute a soluotion in grid points
	assemble_mass_matrixWC6_1D(nx, s, span, r_x, support, K._data)
	# ...
	u   = K.tosparse().dot(control_points)
	# ...
	return u
	
# --------------------------
# In two dimension case
# --------------------------
def CSRBF_basis(X, Y, nx, ny, s):
	'''
	CSRBF_basis : computes the span according to scaling parameter "s" and euclidian distance for all points in the span,
	
	where,
	
	# X  is x-coordinates 2D matrix
	# Y  is y-coordinates 2D matrix
	# nx is number of points in x direction
	# ny is number of points in y direction
	# s  is scaling parameter 
	-> max_span_x : 
	-> max_span_y :
	'''
	support       = zeros((nx,ny), dtype = int)                  # ... Determine the span index for spans
	pyccel_fin_support(X, Y, nx, ny, s, support)
	# ...
	max_span      = int(np.max(support))
	#...
	r_xy          = zeros((nx,ny, max_span))                     # ... computed distance for non niglicted RBF                 
	span          = zeros((nx,ny,4,max_span), dtype = int)    # ... Determine the global index  for non-vanishing CSRBF
	#...
	pyccel_CSRBF_tools(X, Y, nx, ny, s, support, span, r_xy)
	# ...
	max_span_x, max_span_y = np.max(span[:,:,2,:])+1, np.max(span[:,:,3,:])+1
	return max_span_x, max_span_y, span, r_xy, support
	

def results(X, Y, nx, ny, s, control_points, max_span_x = None, max_span_y = None, span = None, r_xy = None, support = None, MQ = None, u_exact = None):
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
	
	# ... Assembles mass matrix for compute a soluotion in grid points
	if MQ is None:
	       #... Mass matrix 
	       K   = StencilMatrix(nx, ny, max_span_x, max_span_y, support, span)
	       assemble_mass_matrixWC6(nx, ny, s, span, r_xy, support, K._data)
	       K = K.tosparse()
	else :
	       #... Mass matrix 
	       K   = zeros((nx,ny,nx,ny), dtype = np.double)
	       assemble_mass_matrixMQ(nx, ny, s, X, Y, K)
	       K   = K.reshape(nx*ny, nx*ny)
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
