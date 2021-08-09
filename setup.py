from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path
import versioneer

here = path.abspath(path.dirname(__file__))
ext_modules = []

# Get the long description from the README file
with open(path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pele_platform",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Automatic platform to launch PELE simulations",
    long_description=long_description,
    url="https://github.com/NostrumBioDiscovery/pele_platform",
    author="Nostrum Biodiscovery",
    author_email="pelesupport@nostrumbiodiscovery.com",
    packages=find_packages(exclude=["docs", "tests", "tests.data"]),
    package_data={
        "pele_platform/AdaptivePELE/atomset": ["*.pxd"],
        "pele_platform/AdaptivePELE/freeEnergies/": ["*.pyx"],
    },
    include_package_data=True,
    install_requires=[
        "PyYAML",
        "fpdf",
        "scikit-learn>0.20",
        "pillow",
        "scipy",
        "matplotlib",
        "mdtraj==1.9.5",
        "biopython",
        "pandas==1.2.3",
        "pytest",
        "cython",
        "numpy",
        "peleffy==1.3.3",
        "hdbscan==0.8.27",
        "AdaptivePELE==1.7.1",
        "frag_pele==3.1.1",
        "PPPele==1.0.8",
        "PlopRotTemp==1.0.1",
        "seaborn"
    ],
    ext_modules=ext_modules,  # accepts a glob pattern
)
