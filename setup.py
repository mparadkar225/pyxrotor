"""Builds the optional compiled hot loops (src/pyxrotor/_speedups.pyx).

The extension is optional: if Cython or a C compiler is missing, or the build fails,
pyxrotor installs without it and runs the pure-Python implementations (identical results,
slower). Floating point must stay IEEE-exact for the results to match XROTOR bit for bit:
no fused multiply-add contraction, no fast-math.
"""
import os

from setuptools import setup
from setuptools.command.build_ext import build_ext

FLAGS = ["-O2", "-ffp-contract=off", "-fno-fast-math"]


class OptionalBuildExt(build_ext):
    def run(self):
        try:
            super().run()
        except Exception as e:                          # noqa: BLE001
            self.warn(f"pyxrotor: compiled speedups not built ({e}); using pure Python")

    def build_extension(self, ext):
        try:
            super().build_extension(ext)
        except Exception as e:                          # noqa: BLE001
            self.warn(f"pyxrotor: compiled speedups not built ({e}); using pure Python")


def extensions():
    if os.environ.get("PYXROTOR_PURE") == "1":
        return []
    try:
        from Cython.Build import cythonize
    except ImportError:
        return []
    from setuptools import Extension
    return cythonize([Extension("pyxrotor._speedups", ["src/pyxrotor/_speedups.pyx"],
                                extra_compile_args=FLAGS, optional=True)],
                     compiler_directives={"language_level": 3})


setup(ext_modules=extensions(), cmdclass={"build_ext": OptionalBuildExt})
