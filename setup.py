from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

extensions = [Extension("meltpack._divergence", ["src/meltpack/_divergence.pyx"],
                        include_dirs=[numpy.get_include()]),
              Extension("meltpack._smooth", ["src/meltpack/_smooth.pyx"],
                        include_dirs=[numpy.get_include()])]

setup(
        name="meltpack",
        version="0.1b1",
        description="Tools for estimating ice shelf melt rates",
        author="Nat Wilson",
        packages=find_packages("src"),
        package_dir={"": "src"},
        ext_modules = cythonize(extensions),
    )
