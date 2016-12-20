from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules=[
    Extension("genereatefracincrements",
              sources=["./cincrements.pyx" ,
                       "./_generatefracincrements.cpp"],
            libraries=["fftw3","m","gsl","gslcblas"], # Unix-like specific
            library_dirs = ['./lib/','/user/lib'],
            include_dirs=[numpy.get_include()],
            language="c++",
	    extra_compile_args=["-O3"]
    )
]
setup(
  name = "cincrements",
  ext_modules = ext_modules,
  cmdclass = {'build_ext': build_ext}
)
