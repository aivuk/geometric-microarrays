import numpy
import ctypes

_libbetag = numpy.ctypeslib.load_library('libbetag', '.')

_libbetag.dist_beta.argtypes = [numpy.ctypeslib.ndpointer(dtype=numpy.float, ndim=2), 
                                    numpy.ctypeslib.ndpointer(dtype=numpy.float, ndim=2),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    numpy.ctypeslib.ndpointer(dtype=numpy.float)]
_libbetag.dist_beta.restype  =  ctypes.c_void_p

def betag(x, y, genes):
    xc = x.ctypes.data_as(numpy.ctypeslib.ndpointer(dtype=numpy.float, ndim=2))
    yc = y.ctypes.data_as(numpy.ctypeslib.ndpointer(dtype=numpy.float, ndim=2)) 
    genesc = genes.ctypes.data_as(numpy.ctypeslib.ndpointer(dtype=numpy.float))  
    nx = x.shape[1]
    ny = y.shape[1] 
    ng = int(len(genes))

    print nx, ny, ng

    _libbetag.dist_beta(xc, yc, 23000, 200, 200, genesc)
    return

x = numpy.random.randn(23000, 200)
y = numpy.random.randn(23000, 200)
genes = numpy.zeros(23000)

betag(x, y, genes)

