import os
from setuptools import setup, find_packages
from setuptools.extension import Extension


def pip_install(pkg_name):
    import subprocess

    subprocess.check_call(
        ["python", "-m", "pip", "install", pkg_name], stdout=subprocess.DEVNULL
    )


try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    pip_install("cython")

    from Cython.Build import cythonize
    from Cython.Distutils import build_ext


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
            "variantpost/pileup.cpp",
            "variantpost/sequence_model.cpp",
            #"variantpost/eval.cpp",
            "variantpost/search.cpp",
            "variantpost/reads.cpp",
            "variantpost/match.cpp",
            #"variantpost/merge.cpp",
            "variantpost/util.cpp",
            #"variantpost/local_reference.cpp",
            #"variantpost/substitutes.cpp",
            #"variantpost/similarity.cpp",
            "variantpost/ssw/ssw.c",
            "variantpost/ssw/ssw_cpp.cpp",
            "variantpost/fasta/Fasta.cpp",
            "variantpost/fasta/split.cpp",
        ],
        language="c++",
        extra_compile_args=["-std=c++17"],
        #extra_compile_args=["-O3", "-std=c++17"],
        extra_link_args=["-std=c++17"],
        #extra_link_args=["-O3", "-std=c++17"],
    ),
]


version = {}
with open("variantpost/version.py") as ver:
    exec(ver.read(), version)

setup(
    name="variantpost",
    version=version["__version__"],
    packages=find_packages(exclude=["tests"]),
    cmdclass={"build_ext": BuildExt},
    ext_modules=cythonize(extensions, annotate=False, language_level="3"),
    install_requires=["pysam>=0.23.3"],
)
