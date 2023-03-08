## coding
#
# Copyright 2023 BAHARI Mustapha
__all__ = ['assemble_Poisson_tools']

from pyccel.decorators              import types
from pyccel.epyccel                 import epyccel


# .... Wendland function C6    
@types('real[:,:]', 'real[:,:]', 'int', 'int', 'real', 'int[:,:,:,:]', 'real[:,:,:]', 'int[:,:]', 'real', 'real[:,:,:,:]', 'real[:,:]')
def assemble_csrbf_tools(X_cor, Y_cor, nx, ny, s, span, r_xy, support, coef_diff, As, Bs): 
     from numpy import exp
     from numpy import cos
     from numpy import sin
     from numpy import pi
     from numpy import sqrt
     
     for i1 in range(nx):
                 # ... on the boundary y = 0.
                 i2 = 0
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1  = span[i1, i2, 0, ij_span]
                         j2  = span[i1, i2, 1, ij_span]
                         r   = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2)
                 # ... Initial conditions
                 x  = X_cor[i1,i2]
                 f = 0.
                 if x <=1.:
                      f = 100.*x
                 else :
                      f = 100.*(2.-x)
                 Bs[i1,i2] = f 
     #Assemble the boundary
     for i2 in range(ny):
                 # ... on the boundary x = 0.
                 i1     = 0
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1  = span[i1, i2, 0, ij_span]
                         j2  = span[i1, i2, 1, ij_span]
                         r   = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2)
                 Bs[i1,i2] = 0.
                 
     for i2 in range(ny):
                 # ... on the boundary x =1
                 i1 = nx-1
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1  = span[i1, i2, 0, ij_span]
                         j2  = span[i1, i2, 1, ij_span]
                         r   = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2)
                 Bs[i1,i2] = 0.


     for i1 in range(1,nx-1):
        for i2 in range(1,ny):
                 
                 x  = X_cor[i1,i2]
                 t  = Y_cor[i1,i2]                 
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1  = span[i1, i2, 0, ij_span]
                         j2  = span[i1, i2, 1, ij_span]
                         r   = r_xy  [i1, i2, ij_span]
                         # ..
                         x2  = X_cor[j1,j2]
                         t2  = Y_cor[j1,j2] 
                         #...
                         As[i1,i2,j1,j2] = -56.*(t-t2)*(1.-r/s)**5/s**2*(1+5.*r/s) - coef_diff * ( (-56/s**2)*(1-r/s)**4.*( (1-r/s)*(1+5*r/s) - 30.*(x-x2)**2/s**2 ) )
                         
                 # ...test 0
                 #f  = 2.*pi**2*sin(pi*x)*sin(pi*y)
                 # ... test 1
                 f = 0.
                 
                 Bs[i1,i2] = f
     return 0
     
assemble_csrbf_stiffnes_rhs = epyccel(assemble_csrbf_tools) 

# ..... MultiQuadric	
@types('real[:,:]', 'real[:,:]', 'int', 'int', 'real', 'int[:,:,:,:]', 'real[:,:,:]', 'int[:,:]', 'real', 'real[:,:,:,:]', 'real[:,:]')
def assemble_rbf_tools(X_cor, Y_cor, nx, ny, c, span, r_xy, support, coef_diff, As, Bs): 
     from numpy import exp
     from numpy import cos
     from numpy import sin
     from numpy import pi
     from numpy import sqrt
     
     for i1 in range(nx):
                 # ... on the boundary y = 0.
                 i2 = 0
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1  = span[i1, i2, 0, ij_span]
                         j2  = span[i1, i2, 1, ij_span]
                         r   = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = sqrt( r**2 + c**2 )
                 # ... Initial conditions
                 x  = X_cor[i1,i2]
                 f = 0.
                 if x <=1.:
                      f = 100.*x
                 else :
                      f = 100.*(2.-x)
                 Bs[i1,i2] = f 
     #Assemble the boundary
     for i2 in range(ny):
                 # ... on the boundary x = 0.
                 i1     = 0
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1  = span[i1, i2, 0, ij_span]
                         j2  = span[i1, i2, 1, ij_span]
                         r   = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = sqrt( r**2 + c**2 )
                 Bs[i1,i2] = 0.
                 
     for i2 in range(ny):
                 # ... on the boundary x =1
                 i1 = nx-1
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1  = span[i1, i2, 0, ij_span]
                         j2  = span[i1, i2, 1, ij_span]
                         r   = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = sqrt( r**2 + c**2 )
                 Bs[i1,i2] = 0.


     for i1 in range(1,nx-1):
        for i2 in range(1,ny):
                 
                 x  = X_cor[i1,i2]
                 t  = Y_cor[i1,i2]                 
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1  = span[i1, i2, 0, ij_span]
                         j2  = span[i1, i2, 1, ij_span]
                         r   = r_xy  [i1, i2, ij_span]
                         # ..
                         x2  = X_cor[j1,j2]
                         t2  = Y_cor[j1,j2] 
                         #...
                         As[i1,i2,j1,j2] = (t-t2)/sqrt( r**2+c**2) - coef_diff * ((t-t2)**2 + c**2)/((sqrt( r**2 + c**2))**3)
                         
                 # ...test 0
                 #f  = 2.*pi**2*sin(pi*x)*sin(pi*y)
                 # ... test 1
                 f = 0.
                 
                 Bs[i1,i2] = f
     return 0
     
assemble_rbf_stiffnes_rhs = epyccel(assemble_rbf_tools) 
