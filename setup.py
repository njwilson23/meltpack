from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

CORR_OBJECTS = ["cross.o", "eval.o", "fft2d.o", "fitreg.o", "gcorr.o",
                "gnorm.o", "kvert.o", "sums.o"]

extensions = [Extension("MELTPACK.correlate", ["src/MELTPACK/correlate.pyx"],
                         extra_objects=CORR_OBJECTS,
                         include_dirs=[numpy.get_include(),
                                       "/usr/lib/gcc/x86_64-linux-gnu/5"],
                         libraries=["gfortran"]),
               Extension("MELTPACK._divergence", ["src/MELTPACK/_divergence.pyx"],
                          include_dirs=[numpy.get_include()])
]

setup(
        name="MELTPACK",
        version="0.1a2",
        description="Tools for estimating ice shelf melt rates",
        author="Nat Wilson",
        packages=find_packages("src"),
        package_dir={"": "src"},

        # ext_modules = cythonize("src/MELTPACK/*.pyx"),
        ext_modules = cythonize(extensions),
    )
