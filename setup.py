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
    Extension(
        "variantpost.preprocessor", ["variantpost/preprocessor.pyx"], language="c++",
    ),
    Extension(
        "variantpost.processor_wrapper",
        [
            "variantpost/processor_wrapper.pyx",
            "variantpost/processor.cpp",
            "variantpost/pileup_parser.cpp",
            "variantpost/util.cpp",
            "variantpost/swlib.cpp",
            "variantpost/localn.cpp",
            "variantpost/ssw/ssw.c",
            "variantpost/ssw/ssw_cpp.cpp",
            "variantpost/fasta/Fasta.cpp",
            "variantpost/fasta/split.cpp",
            "variantpost/aligned_variant.cpp",
        ],
        language="c++",
    ),
]

setup(cmdclass={"build_ext": BuildExt}, ext_modules=cythonize(extensions, annotate=False))
