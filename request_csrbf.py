## coding
#
# Copyright 2023 BAHARI Mustapha

"""
Basic module that provides the means for evaluating the CSRBF basis
functions. In order to simplify automatic Fortran code
generation with Pyccel.
"""
__all__ = ['pyccel_fin_support_1D',
          'pyccel_CSRBF_tools_1D',
          'assemble_mass_matrixWC6_1D',
          'pyccel_fin_support',
          'pyccel_CSRBF_tools',
          'assemble_mass_matrixWC6',
          'assemble_mass_matrixMQ',]

from pyccel.decorators              import types
from pyccel.epyccel                 import epyccel

# --------------------------
# In one dimension case
# --------------------------
@types('real[:]', 'int', 'real', 'int[:]')
def fin_support_1D(X, nx, s, support):

	for i in range(nx):
	        i_span = 0
	        for k1 in range(1,i+1):
	            r = abs(X[i]-X[i-k1])
	            if r < s :
	                 i_span            += 1
	        for k1 in range(0,nx-i):
	            r = abs(X[i]-X[i+k1])
	            if r < s :
	                 i_span            += 1
	        support[i] = i_span
	return 0.
pyccel_fin_support_1D = epyccel(fin_support_1D)

# -------
@types('real[:]', 'int', 'real', 'int[:]', 'int[:,:]', 'real[:,:]')
def CSRBF_tools_1D(X, nx, s, support, span, r_span):

	for i in range(nx):

	        i_span   = 0
	        max_span = support[i]
	        for k1 in range(0,nx-i):
	            r = abs(X[i]-X[i+k1])
	            if r < s :
	                 span[i,i_span]   = i+k1

	                 r_span[i,i_span] = r
	                 
	                 i_span          += 1
	            if i_span == max_span:
	                    break
	        for k1 in range(1,i+1):
	            r = abs(X[i]-X[i-k1])
	            if r < s :
	                 span[i,i_span]   = i-k1
	                 
	                 r_span[i,i_span] = r
	                 i_span          += 1
	            if i_span == max_span:
	                    break
	return 0.
# ...
pyccel_CSRBF_tools_1D = epyccel(CSRBF_tools_1D)

	
# .... Wendland function C6
@types('int', 'real', 'int[:,:]', 'real[:,:]', 'int[:]', 'real[:,:]')
def assemble_mass_1D(nx, s, span, r_x, support, K): 

     for i1 in range(0,nx):
                 
            spectr = support[i1]
            for i_span in range(spectr):
                r  = r_x[i1, i_span]
                #...
                K[i1, i_span] = (1-r/s)**6.*(3+18.*r/s+35.*(r/s)**2)
     return 0
     
assemble_mass_matrixWC6_1D = epyccel(assemble_mass_1D)    

# --------------------------
# In two dimension case
# --------------------------
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
@types('real[:,:]', 'real[:,:]', 'int', 'int', 'real', 'int[:,:]', 'int[:,:,:,:]', 'real[:,:,:]')
def CSRBF_tools(X, Y, nx, ny, s, support, span, r_span):
	from numpy import sqrt
	for i in range(nx):
	  for j in range(ny):
	    
	    ij_span  = 0
	    max_span = support[i,j]
	    i_loc    = 0
	    for k1 in range(0,nx-i):
	        j_loc  = 0
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i+k1,j+k2])**2 + (Y[i,j]-Y[i+k1,j+k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span] = i+k1
	                 span[i,j,1,ij_span] = j+k2 	
	                                  
	                 span[i,j,2,ij_span] = i_loc	                 
	                 span[i,j,3,ij_span] = j_loc
	                 
	                 r_span[i,j,ij_span] = r
	                 
	                 ij_span            += 1
	                 j_loc              += 1
	            if ij_span == max_span:
	                    break
	        for k2 in range(1,j+1):
	            r = sqrt((X[i,j]-X[i+k1,j-k2])**2 + (Y[i,j]-Y[i+k1,j-k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span] = i+k1
	                 span[i,j,1,ij_span] = j-k2

	                 span[i,j,2,ij_span] = i_loc	                 
	                 span[i,j,3,ij_span] = j_loc
	                 
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1
	                 j_loc              += 1
	            if ij_span == max_span:
	                    break
	        if j_loc > 0 :
	               i_loc += 1
	        if ij_span == max_span:
	               break

	    for k1 in range(1,i+1):
	        j_loc  = 0
	        for k2 in range(0,ny-j):
	            r = sqrt((X[i,j]-X[i-k1,j+k2])**2 + (Y[i,j]-Y[i-k1,j+k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span] = i-k1
	                 span[i,j,1,ij_span] = j+k2

	                 span[i,j,2,ij_span] = i_loc	                 
	                 span[i,j,3,ij_span] = j_loc
	                 
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1
	                 j_loc              += 1
	            if ij_span == max_span:
	                    break
	        for k2 in range(1,j+1):
	            r = sqrt((X[i,j]-X[i-k1,j-k2])**2 + (Y[i,j]-Y[i-k1,j-k2])**2)
	            if r < s :
	                 span[i,j,0,ij_span] = i-k1
	                 span[i,j,1,ij_span] = j-k2

	                 span[i,j,2,ij_span] = i_loc	                 
	                 span[i,j,3,ij_span] = j_loc
	                 
	                 r_span[i,j,ij_span] = r
	                 ij_span            += 1
	                 j_loc              += 1
	            if ij_span == max_span:
	                    break
	        if j_loc > 0 :
	               i_loc += 1
	        if ij_span == max_span:
	               break
	return 0.
# ...
pyccel_CSRBF_tools = epyccel(CSRBF_tools)

	
# .... Wendland function C6
@types('int', 'int', 'real', 'int[:,:,:,:]', 'real[:,:,:]', 'int[:,:]', 'real[:,:,:,:]',)
def assemble_mass(nx, ny, s, span, r_xy, support, K): 

     for i1 in range(0,nx):
        for i2 in range(0,ny):
                 
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1 = span[i1, i2, 2, ij_span]
                         j2 = span[i1, i2, 3, ij_span]
                         r  = r_xy  [i1, i2, ij_span]
                         #...
                         K[i1,i2,j1,j2]  = (1-r/s)**6.*(3+18.*r/s+35.*(r/s)**2)
     return 0
     
assemble_mass_matrixWC6 = epyccel(assemble_mass)    


