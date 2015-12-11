from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
import numpy

OBJECTS = ["melt.o", "w2f__types.o", "OAD_active.o", "OAD_cp.o", "OAD_tape.o",
           "OAD_rev.o", "melt_adjoint.o", "adjoint_driver.o",
           "utility.o", "cwrap_melt.o"]

setup(
        name="MELTPACK",
        version="0.1a",
        description="Tools for estimating ice shelf melt rates",
        author="Nat Wilson",
        packages=find_packages("src"),
        package_dir={"": "src"},

        #cmdclass = {'build_ext': build_ext},
        #ext_modules = [Extension("meltfuncs",
        #                         sources=["meltfuncs.pyx"],
        #                         extra_objects=OBJECTS,
        #                         extra_compile_args=["-fopenmp", "-Wall"],
        #                         extra_link_args=["-fopenmp"],
        #                         include_dirs=[numpy.get_include(),
        #                                       "/usr/lib/gcc/x86_64-linux-gnu/5"],
        #                         libraries=["gfortran"])],
    )
