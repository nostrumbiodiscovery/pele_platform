import pele_platform
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from distutils.extension import Extension

here = path.abspath(path.dirname(__file__))
ext_modules = []
cmdclass = {}
# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
     long_description = f.read()

setup(
    name="pele_platform",
    version=pele_platform.__version__,
    description='Automatic platform to launch PELE',
    long_description=long_description,
    url="https://github.com/NostrumBioDiscovery/pele_platform",
    author='Daniel Soler',
    author_email='daniel.soler@nostrumbiodiscovery.com',
    packages=find_packages(exclude=['docs', 'tests']),
    package_data={"pele_platform/AdaptivePELE/atomset": ['*.pxd'], "pele_platform/AdaptivePELE/freeEnergies/": ['*.pyx']},
    include_package_data=True,
    install_requires=["AdaptivePELE", "PPPele", "PyYAML", "fpdf", "scikit-learn", "pillow", "scipy", "matplotlib", 
       "biopython", "pandas", "pytest", "cython", "numpy"],    
    cmdclass=cmdclass,
    ext_modules=ext_modules  # accepts a glob pattern
)
