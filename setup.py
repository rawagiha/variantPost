from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import setup, find_packages
from setuptools.extension import Extension


class BuildExt(build_ext):
    def build_extensions(self):
        if "-Wstrict-prototypes" in self.compiler.compiler_so:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        super().build_extensions()


extensions = [
    Extension(
        "variantpost.__search",
        [
            "variantpost/__search.pyx",
            "variantpost/contig.cpp",
            "variantpost/eval.cpp",
            "variantpost/search.cpp",
            "variantpost/reads.cpp",
            "variantpost/match.cpp",
            "variantpost/merge.cpp",
            "variantpost/util.cpp",
            "variantpost/substitutes.cpp",
            "variantpost/ssw/ssw.c",
            "variantpost/ssw/ssw_cpp.cpp",
            "variantpost/fasta/Fasta.cpp",
            "variantpost/fasta/split.cpp",
        ],
        language="c++",
        extra_compile_args=["-O3", "-std=c++17"],
        extra_link_args=["-O3", "-std=c++17"],
    ),
]

setup(
    name="variantpost",
    packages=find_packages(exclude=["tests"]),
    cmdclass={"build_ext": BuildExt},
    ext_modules=cythonize(extensions, annotate=False, language_level="3"),
)
