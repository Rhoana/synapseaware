from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        name='wiring',
        include_dirs=[np.get_include()],
        sources=['wiring.pyx', 'cpp-wiring.cpp', 'cpp-thinning.cpp', 'cpp-refinement.cpp', 'cpp-MinBinaryHeap.cpp'],
        extra_compile_args=['-O4', '-std=c++0x'],
        language='c++'
    )
]

setup(
    name='wiring',
    ext_modules=cythonize(extensions)
)
