from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        name='topological_thinning',
        include_dirs=[np.get_include()],
        sources=['topological_thinning.pyx', 'cpp-topological-thinning.cpp'],
        extra_compile_args=['-O4', '-std=c++0x'],
        language='c++'
    )
]

setup(
    name='topological_thinning',
    ext_modules=cythonize(extensions)
)
