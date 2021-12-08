from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import setup
from setuptools.extension import Extension


class BuildExt(build_ext):
    def build_extensions(self):
        if "-Wstrict-prototypes" in self.compiler.compiler_so:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        super().build_extensions()


extensions = [
    Extension(
        "variantpost.variantalignment",
        ["variantpost/variantalignment.pyx"],
        language="c++",
    ),
    #Extension("variantpost.util", ["variantpost/util.pyx"], language="c++",),

]

setup(cmdclass={"build_ext": BuildExt}, ext_modules=cythonize(extensions))
