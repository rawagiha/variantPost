from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import setup, find_packages
from setuptools.extension import Extension
#from pysam import get_include as pysam_get_include

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
        #include_dirs=pysam_get_include(),
        extra_compile_args=["-O3", "-std=c++17"],
        extra_link_args=["-O3", "-std=c++17"],
    ),
    #Extension(
    #    "variantpost.preprocessor", ["variantpost/preprocessor.pyx"], language="c++",
    #    extra_compile_args=["-O3", "-std=c++17"],
    #    extra_link_args=["-O3", "-std=c++17"],
    #),
    Extension(
        "variantpost.cy_search",
        [
            "variantpost/cy_search.pyx",
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
        #include_dirs=pysam_get_include(),
        extra_compile_args=["-O3", "-std=c++17"],
        extra_link_args=["-O3", "-std=c++17"],
    ),
]

setup(name="variantpost",  packages=find_packages(exclude=["tests"]), cmdclass={"build_ext": BuildExt}, ext_modules=cythonize(extensions, annotate=False, language_level="3"),)
