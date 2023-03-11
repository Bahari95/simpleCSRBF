## coding
#
# Copyright 2023 BAHARI Mustapha

"""
assemble_sol_exact : code that compute the analytical soluion of Fourier equation accelerated with Pyccel.
assemble_massMQ    : code assembles multiquadratic RBF
"""
__all__ = ['exact_sol',
          'results']
from pyccel.decorators              import types
from pyccel.epyccel                 import epyccel

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
exact_sol = epyccel(assemble_sol_exact) 


# ..... MultiQuadric in 2D
@types('int', 'int', 'real', 'real[:,:]', 'real[:,:]', 'real[:,:,:,:]')
def assemble_massMQ(nx, ny, c, X_cor, Y_cor, K): 

     from numpy import sqrt
     for i1 in range(0,nx):
        for i2 in range(0,ny):
                 
           for j1 in range(0,nx):
                 for j2 in range(0,ny):
                         r  = (X_cor[i1, i2]-X_cor[j1, j2])**2+(Y_cor[i1, i2]-Y_cor[j1, j2])**2
                         #...
                         K[i1,i2,j1,j2]  = sqrt(r + c**2)
     return 0
     
assemble_mass_matrixMQ = epyccel(assemble_massMQ)    
