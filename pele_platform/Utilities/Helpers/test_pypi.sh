#Clean after build
pip uninstall -r requirements.txt --yes
pip uninstall pele_platform --yes
pip install numpy cython
pip install pele_platform
cd tests
pytest


