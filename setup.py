import os
import sys
from setuptools import setup, find_packages
from setuptools.extension import Extension

version = {}
with open("variantpost/version.py", encoding="utf-8") as ver:
    exec(ver.read(), version)

try:
    from Cython.Distutils import build_ext as _build_ext
except ImportError:
    from setuptools.command.build_ext import build_ext as _build_ext


class SafeBuildExt(_build_ext):
    def build_extensions(self):
        if hasattr(self.compiler, "compiler_so") and isinstance(
            self.compiler.compiler_so, list
        ):
            if "-Wstrict-prototypes" in self.compiler.compiler_so:
                self.compiler.compiler_so.remove("-Wstrict-prototypes")

        import pysam

        pysam_includes = pysam.get_include()
        if isinstance(pysam_includes, str):
            pysam_includes = [pysam_includes]
        htslib_includes = [os.path.join(p, "htslib") for p in pysam_includes]

        for ext in self.extensions:
            ext.include_dirs.extend(pysam_includes + htslib_includes)

        super().build_extensions()


extra_compile_args = []
if sys.platform == "win32":
    extra_compile_args = ["/std:c++17"]
else:
    extra_compile_args = ["-std=c++17"]


extensions = [
    Extension(
        "variantpost.__search",
        [
            "variantpost/__search.pyx",
            "variantpost/pileup.cpp",
            "variantpost/search.cpp",
            "variantpost/reads.cpp",
            "variantpost/match.cpp",
            "variantpost/util.cpp",
            "variantpost/consensus.cpp",
            "variantpost/ssw/ssw.c",
            "variantpost/ssw/ssw_cpp.cpp",
            "variantpost/fasta/Fasta.cpp",
            "variantpost/fasta/split.cpp",
        ],
        language="c++",
        include_dirs=[],
        extra_compile_args=extra_compile_args,
    ),
]


def get_ext_modules():
    try:
        from Cython.Build import cythonize

        return cythonize(extensions, annotate=False, language_level="3")
    except ImportError:
        return extensions


setup(
    name="variantpost",
    version=version["__version__"],
    packages=find_packages(exclude=["tests"]),
    cmdclass={"build_ext": SafeBuildExt},
    ext_modules=get_ext_modules(),
    setup_requires=["pysam>=0.23.3", "cython>=3.0.0"],
    install_requires=["pysam>=0.23.3"],
    entry_points={
        "console_scripts": [
            "indelinside=variantpost.__main__:main",
        ],
    },
)
