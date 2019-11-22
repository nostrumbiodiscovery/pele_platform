#Clean and reinstall requirements
pip uninstall -r requirements.txt --yes
pip uninstall mdtraj --yes
pip uninstall pele_platform --yes
pip install cython numpy
pip install .

#fast test
python -m pele_platform.main -h

#Clean and build
rm -r dist pele_platform.egg*
python setup.py sdist
twine upload  --repository-url https://test.pypi.org/legacy/ dist/*
#twine upload  dist/*

