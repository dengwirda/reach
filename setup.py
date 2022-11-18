
import io
import os
from setuptools import setup, find_packages
from setuptools.extension import Extension

# python3 setup.py build_ext --inplace

# https://stackoverflow.com/questions/4505747/
#   how-should-i-structure-a-python-package-that-contains-cython-code

# https://stackoverflow.com/questions/14657375/
#   cython-fatal-error-numpy-arrayobject-h-no-such-file-or-directory

# https://github.com/cython/cython/issues/1480

# https://stackoverflow.com/questions/14372706/
#   visual-studio-cant-build-due-to-rc-exe

EXT_MODULES = []

try:
    from Cython.Build import cythonize
    EXT_MODULES += cythonize(Extension(
        "reach.kernel",
        sources=[os.path.join("reach", "kernel.pyx")],
        include_dirs=[])
    )

except ImportError:
    EXT_MODULES += [Extension(
        "reach.kernel",
        sources=[os.path.join("reach", "kernel.c")],
        include_dirs=[])
    ]

NAME = "reach"
DESCRIPTION = "River network simplification."
AUTHOR = "Darren Engwirda"
AUTHOR_EMAIL = "dengwirda@lanl.gov"
URL = "https://github.com/dengwirda/reach"
VERSION = "0.0.1"
REQUIRES_PYTHON = ">=3.3.0"
KEYWORDS = "Hydrology Earth System Modelling"

REQUIRED = [
    "numpy", "netCDF4", "argparse"
]

CLASSIFY = [
    "Development Status :: 4 - Beta",
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Mathematics",
    "Topic :: Scientific/Engineering :: Physics",
    "Topic :: Scientific/Engineering :: GIS"
]

HERE = os.path.abspath(os.path.dirname(__file__))

try:
    with io.open(os.path.join(
            HERE, "README.md"), encoding="utf-8") as f:
        LONG_DESCRIPTION = "\n" + f.read()

except FileNotFoundError:
    LONG_DESCRIPTION = DESCRIPTION

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    license="custom",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    python_requires=REQUIRES_PYTHON,
    keywords=KEYWORDS,
    url=URL,
    packages=find_packages(),
    ext_modules=EXT_MODULES,
    install_requires=REQUIRED,
    classifiers=CLASSIFY
)
