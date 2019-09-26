import platform
from distutils.core import setup
from distutils.extension import Extension

import numpy
from Cython.Build import cythonize
from Cython.Distutils import build_ext

compile_extra_args = []
link_extra_args = []

if platform.system() == "Darwin":
    compile_extra_args = ['-std=c++11', "-mmacosx-version-min=10.9"]
    link_extra_args = ["-stdlib=libc++", "-mmacosx-version-min=10.9"]

ext_modules = [Extension("inclose5_path_cython", sources=["inclose5_path_cython.pyx"],
                         language='c++',
                         extra_compile_args=compile_extra_args,
                         extra_link_args=link_extra_args,
                         include_dirs=[numpy.get_include()]),
               Extension("set_utils", sources=["set_utils.pyx"],
                         language='c++',
                         extra_compile_args=compile_extra_args,
                         extra_link_args=link_extra_args,
                         include_dirs=[numpy.get_include()])
               ]

setup(ext_modules = cythonize(
        ext_modules,
        compiler_directives={'language_level': "3"},
        annotate=True
))