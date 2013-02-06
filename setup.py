from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

source_files = ['beta.pyx']

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("bg", source_files)]
)

