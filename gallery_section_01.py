## coding
#
# Copyright 2023 BAHARI Mustapha
__all__ = ['assemble_Poisson_tools_1D',
           'assemble_Poisson_tools']

from pyccel.decorators              import types
from pyccel.epyccel                 import epyccel

# --------------------------
# In one dimension case
# --------------------------
@types('real[:]', 'int', 'real', 'int[:,:]', 'real[:,:]', 'int[:]', 'real[:,:]', 'real[:]')
def assemble_Poisson_tools_1D(X_cor, nx, s, span, r_x, support, A, B): 
     from numpy import exp
     from numpy import cos
     from numpy import sin
     from numpy import pi
     from numpy import sqrt

     for i1 in range(1,nx-1):
                 
          spectr = support[i1]
          for i_span in range(spectr):
              r  = r_x[i1, i_span]
              #...
              A[i1,i_span] = -(-56/s**2)*(1-r/s)**4 *(1+4.*r/s-35.*(r/s)**2)
                         
          x  = X_cor[i1]
          # ...test 0
          f  = pi**2*sin(pi*x)
                 
          B[i1] = f

     # ... on the boundary x = 0.
     i1     = 0
     spectr = support[i1]
     for i_span in range(spectr):

             r  = r_x[i1, i_span]
             #...
             A[i1,i_span] = (1-r/s)**6 *(3+18*r/s+35.*(r/s)**2)
     B[i1] = 0
                 
     # ... on the boundary x =1
     i1     = nx-1
     spectr = support[i1]
     for i_span in range(spectr):
          
          #..
          r  = r_x[i1, i_span]
          #...
          A[i1,i_span] = (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2)
     B[i1] = 0
     return 0
     
assemble_Poisson_stiffnes_rhs_1D = epyccel(assemble_Poisson_tools_1D) 

# --------------------------
# In two dimension case
# --------------------------
@types('real[:,:]', 'real[:,:]', 'int', 'int', 'real', 'int[:,:,:,:]', 'real[:,:,:]', 'int[:,:]', 'real[:,:,:,:]', 'real[:,:]')
def assemble_Poisson_tools(X_cor, Y_cor, nx, ny, s, span, r_xy, support, As, Bs): 
     from numpy import exp
     from numpy import cos
     from numpy import sin
     from numpy import pi
     from numpy import sqrt

     for i1 in range(1,nx-1):
        for i2 in range(1,ny-1):
                 
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         j1 = span[i1, i2, 2, ij_span]
                         j2 = span[i1, i2, 3, ij_span]
                         r  = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = -(-56/s**2)*(1-r/s)**4.*(2+8.*r/s-40.*(r/s)**2)
                         
                 x  = X_cor[i1,i2]
                 y  = Y_cor[i1,i2]
                 # ...test 0
                 #f  = -(15625.0 - 31250.0*x)*(2*x - 1.0)*exp(-7812.5*((x - 0.5)**2 + (y - 0.5)**2 - 0.096)**2) 
                 #f += 7812.5*(15625.0 - 31250.0*x)*(4*x - 2.0)*((x - 0.5)**2 + (y - 0.5)**2 - 0.096)**2*exp(-7812.5*((x - 0.5)**2 + (y - 0.5)**2 - 0.096)**2) 
                 #f += -(15625.0 - 31250.0*y)*(2*y - 1.0)*exp(-7812.5*((x - 0.5)**2 + (y - 0.5)**2 - 0.096)**2) 
                 #f += 7812.5*(15625.0 - 31250.0*y)*(4*y - 2.0)*((x - 0.5)**2 + (y - 0.5)**2 - 0.096)**2*exp(-7812.5*((x - 0.5)**2 + (y - 0.5)**2 - 0.096)**2) 
                 #f += -2*(-31250.0*(x - 0.5)**2 - 31250.0*(y - 0.5)**2 + 3000.0)*exp(-7812.5*((x - 0.5)**2 + (y - 0.5)**2 - 0.096)**2)
                 # ...test 1
                 f  = 2.*pi**2*sin(pi*x)*sin(pi*y)
                 
                 Bs[i1,i2] = f

     #Assemble the boundary
     for i2 in range(ny):
                 # ... on the boundary x = 0.
                 i1     = 0
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):

                         j1 = span[i1, i2, 2, ij_span]
                         j2 = span[i1, i2, 3, ij_span]
                         r  = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2)
                 Bs[i1,i2] = 0
                 
     for i2 in range(ny):
                 # ... on the boundary x =1
                 i1 = nx-1
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         #..
                         j1 = span[i1, i2, 2, ij_span]
                         j2 = span[i1, i2, 3, ij_span]
                         r  = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2)
                 Bs[i1,i2] = 0
     for i1 in range(nx):
                 # ... on the boundary y = 0.
                 i2 = 0
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         #...
                         j1 = span[i1, i2, 2, ij_span]
                         j2 = span[i1, i2, 3, ij_span]
                         r  = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2)
                 Bs[i1,i2] = 0
                 
     for i1 in range(nx):
                 # ... on the boundary y =1
                 i2 = ny-1
                 spectr = support[i1, i2]
                 for ij_span in range(spectr):
                         #...
                         j1 = span[i1, i2, 2, ij_span]
                         j2 = span[i1, i2, 3, ij_span]
                         r  = r_xy  [i1, i2, ij_span]
                         #...
                         As[i1,i2,j1,j2] = (1-r/s)**6.*(3+18*r/s+35.*(r/s)**2)
                 Bs[i1,i2] = 0
     return 0
     
assemble_Poisson_stiffnes_rhs = epyccel(assemble_Poisson_tools) 
