"""Builds pyxrotor's compiled core, src/pyxrotor/_speedups.pyx (Cython; needs a C compiler).

Floating point must stay IEEE-exact for the results to match XROTOR bit for bit: no fused
multiply-add contraction, no fast-math.
"""
from Cython.Build import cythonize
from setuptools import Extension, setup

FLAGS = ["-O2", "-ffp-contract=off", "-fno-fast-math"]

setup(ext_modules=cythonize(
    [Extension("pyxrotor._speedups", ["src/pyxrotor/_speedups.pyx"], extra_compile_args=FLAGS)],
    compiler_directives={"language_level": 3}))
