
import numpy
from cffi import FFI
import ctypes

_libbetag = numpy.ctypeslib.load_library('libbetag', '.')

_libbetag.dist_beta.argtypes = [numpy.ctypeslib.ndpointer(dtype=numpy.float, ndim=2), 
                                    numpy.ctypeslib.ndpointer(dtype=numpy.float, ndim=2),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    numpy.ctypeslib.ndpointer(dtype=numpy.float,
                                        ndim=1)]
_libbetag.dist_beta.restype  =  ctypes.c_void_p

ffi = FFI()

def betag(x, y, genes):
    xc = numpy.asarray(x, dtype=numpy.float)
    yc = numpy.asarray(y, dtype=numpy.float)
    genesc = numpy.asarray(genes, dtype=numpy.float)
    _libbetag.dist_beta(xc, yc, int(len(genes)), int(x.shape[1]), int(y.shape[1]), genesc)
    return

x = numpy.random.randn(23000, 200)
y = numpy.random.randn(23000, 200)
genes = numpy.zeros(23000)

betag(x, y, genes)
