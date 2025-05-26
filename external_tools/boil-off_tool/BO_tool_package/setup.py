from distutils.core import setup
from Cython.Build import cythonize
import numpy


setup(ext_modules = cythonize(r'C:\Users\hp\Desktop\University\DSE\DSE-07\boil-off_tool\BO_tool_package\BO_tool_main_body.pyx'), include_dirs=[numpy.get_include()])

## Remember to run the line below in the Terminal every time that changes are made to the *.pyx file##
##### python setup.py build_ext --inplace #####
## This is to convert the *.pyx file into *.c and *so, so that Cython can run it##