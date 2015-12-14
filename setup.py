from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

setup(
        name="MELTPACK",
        version="0.1a1",
        description="Tools for estimating ice shelf melt rates",
        author="Nat Wilson",
        packages=find_packages("src"),
        package_dir={"": "src"},

        ext_modules = cythonize("src/MELTPACK/*.pyx"),
    )
