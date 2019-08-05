from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        name='teaser',
        include_dirs=[np.get_include()],
        sources=['teaser.pyx', 'cpp-teaser.cpp', 'cpp-MinBinaryHeap.cpp'],
        extra_compile_args=['-O4', '-std=c++0x'],
        language='c++'
    )
]

setup(
    name='teaser',
    ext_modules=cythonize(extensions)
)
