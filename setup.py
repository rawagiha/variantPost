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
    Extension("variantpost.preprocessor", ["variantpost/preprocessor.pyx"], language="c++",),
    Extension("variantpost.pileup_parser_wrapper", ["variantpost/pileup_parser_wrapper.pyx", "variantpost/pileup_parser.cpp", "variantpost/util.cpp", "variantpost/swlib.cpp", "variantpost/ssw/ssw.c", "variantpost/ssw/ssw_cpp.cpp"], language="c++",),

]

setup(cmdclass={"build_ext": BuildExt}, ext_modules=cythonize(extensions))
