#
# Copyright 2023 BAHARI MUSTAPHA

import numpy as np
from scipy.sparse import coo_matrix

#===============================================================================
class StencilMatrix( object ):
    """
    Matrix in n-dimensional stencil format.
    """
    def __init__( self, nx, ny, max_span_x, max_span_y, support, span, data_type = None):

        #assert isinstance( dim, 2 )
        
        diags           = (max_span_x, max_span_y)
        dims            = (nx, ny)
        self._support   = support
        self._ndim      = dims
        self._span      = span
        
        if data_type is None:
           data_type    = np.float64
           
        self._data_type = data_type
        self._data      = np.zeros( dims+diags, dtype= data_type )
        
    # ...
    def tosparse( self):

        coo = self._tocoo_no_pads()

        return coo
    #...
    def _tocoo_no_pads( self ):

        # Shortcuts
        npoints = self._support
        nx, ny  = self._ndim
        span    = self._span
        dtype   = self._data_type
        # COO storage
        rows = []
        cols = []
        data = []
        for i1 in range(0,nx):
            for i2 in range(0,ny):
                 
                 spectr = npoints[i1, i2]
                 for ij_span in range(spectr):
                         j1    = span[i1, i2, 2, ij_span]
                         j2    = span[i1, i2, 3, ij_span]
                         # ...
                         #print(i1,i2,j1,j2)
                         value = self._data[i1, i2, j1, j2]
                         # ...
                         j1    = span[i1, i2, 0, ij_span]
                         j2    = span[i1, i2, 1, ij_span]
                         I     = i1*ny + i2
                         J     = j1*ny + j2 
                         rows.append( I )
                         cols.append( J )
                         data.append( value)
        M = coo_matrix(
                (data,(rows,cols)),
                shape = [nx*ny,nx*ny],
                dtype = dtype
        )

        M.eliminate_zeros()
        return M
