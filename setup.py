import os
import pysam
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

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

pysam_includes = pysam.get_include()
if isinstance(pysam_includes, str):
    pysam_includes = [pysam_includes]

htslib_includes = [os.path.join(p, "htslib") for p in pysam_includes]

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
        include_dirs=pysam_includes + htslib_includes,
        extra_compile_args=["-std=c++17"],
        extra_link_args=["-std=c++17"],
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

    entry_points={
         'console_scripts': [
             'indelinside=variantpost.__main__:main',
         ],
     },
)
