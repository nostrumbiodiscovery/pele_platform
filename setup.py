import numpy
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from distutils.extension import Extension
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True
    print(use_cython)
from distutils.command.sdist import sdist as _sdist

# Run the following line to compile atomset package
# python setup.py build_ext --inplace


class sdist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are
        # up-to-date
        from Cython.Build import cythonize
        cythonize(['cython/mycythonmodule.pyx'])
        _sdist.run(self)
        cmdclass['sdist'] = sdist

here = path.abspath(path.dirname(__file__))
ext_modules = []
cmdclass = {}
# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
     long_description = f.read()

if use_cython:
    ext_modules += [
            Extension("pele_platform.AdaptivePELE.atomset.atomset", ["pele_platform/AdaptivePELE/atomset/atomset.pyx"], include_dirs=["pele_platform/AdaptivePELE", "pele_platform/AdaptivePELE/atomset"]),
        Extension("pele_platform.AdaptivePELE.atomset.SymmetryContactMapEvaluator", ["pele_platform/AdaptivePELE/atomset/SymmetryContactMapEvaluator.pyx"], include_dirs=["pele_platform/AdaptivePELE", "pele_platform/AdaptivePELE/atomset"]),
        Extension("pele_platform.AdaptivePELE.atomset.RMSDCalculator", ["pele_platform/AdaptivePELE/atomset/RMSDCalculator.pyx"], include_dirs=["pele_platform/AdaptivePELE", "pele_platform/AdaptivePELE/atomset"]),
        Extension("pele_platform.AdaptivePELE.freeEnergies.utils", ["pele_platform/AdaptivePELE/freeEnergies/utils.pyx"], include_dirs=["pele_platform/AdaptivePELE", "pele_platform/AdaptivePELE/freeEnergies"])
    ]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("pele_platform.AdaptivePELE.atomset.atomset", ["pele_platform/AdaptivePELE/atomset/atomset.c"], include_dirs=["pele_platform/AdaptivePELE", "pele_platform/AdaptivePELE/atomset"]),
        Extension("pele_platform.AdaptivePELE.atomset.SymmetryContactMapEvaluator", ["pele_platform/AdaptivePELE/atomset/SymmetryContactMapEvaluator.c"], include_dirs=["pele_platform/AdaptivePELE", "pele_platform/Utilities/AdaptivePELE/atomset"]),
        Extension("pele_platform.AdaptivePELE.atomset.RMSDCalculator", ["pele_platform/AdaptivePELE/atomset/RMSDCalculator.c"], include_dirs=["pele_platform/AdaptivePELE", "pele_platform/AdaptivePELE/atomset"]),
        Extension("pele_platform.AdaptivePELE.freeEnergies.utils", ["pele_platform/AdaptivePELE/freeEnergies/utils.c"], include_dirs=["pele_platform/AdaptivePELE", "pele_platform/AdaptivePELE/freeEnergies"])
    ]

setup(
    name="pele_platform",
    version="1.0.0.1",
    description='Automatic platform to launch PELE',
    long_description=long_description,
    url="https://github.com/NostrumBioDiscovery/pele_platform",
    author='Daniel Soler',
    author_email='daniel.soler@nostrumbiodiscovery.com',
    packages=find_packages(exclude=['docs', 'tests']),
    package_data={"pele_platform/AdaptivePELE/atomset": ['*.pxd']},
    include_package_data=True,
    include_dirs=[numpy.get_include()],
    install_requires=['cython', 'numpy', 'scipy', 'matplotlib', 'biopython ', 'pandas', 'pyemma', 'prody', 'six', 'future', 'fpdf', 'pytest', 'mdtraj'],
    cmdclass=cmdclass,
    ext_modules=cythonize(ext_modules)  # accepts a glob pattern
)
