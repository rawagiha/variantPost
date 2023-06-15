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
        extra_compile_args=['-O3'],
        extra_link_args=['-O3']
    ),
    Extension(
        "variantpost.preprocessor", ["variantpost/preprocessor.pyx"], language="c++",
        extra_compile_args=['-O3'],
        extra_link_args=['-O3']
    ),
    Extension(
        "variantpost._wrapper",
        [
            "variantpost/_wrapper.pyx",
            #"variantpost/processor.cpp",
            "variantpost/eval.cpp",
            "variantpost/search.cpp",
            #"variantpost/pileup_parser.cpp",
            "variantpost/reads.cpp",
            "variantpost/merge.cpp",
            "variantpost/match.cpp",
            #"variantpost/read_classifier.cpp",
            "variantpost/util.cpp",
            "variantpost/contig.cpp",
            #"variantpost/swlib.cpp",
            #"variantpost/localn.cpp",
            "variantpost/ssw/ssw.c",
            "variantpost/ssw/ssw_cpp.cpp",
            "variantpost/fasta/Fasta.cpp",
            "variantpost/fasta/split.cpp",
            #"variantpost/aligned_target.cpp",
            #"variantpost/unaligned_target.cpp",
        ],
        language="c++",
        extra_compile_args=['-O3'],
        extra_link_args=['-O3']
    ),
]

setup(cmdclass={"build_ext": BuildExt}, ext_modules=cythonize(extensions, annotate=False))
