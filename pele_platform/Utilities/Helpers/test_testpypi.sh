#Clean after build
pip uninstall -r requirements.txt --yes
pip uninstall pele_platform --yes
pip install numpy cython
pip install -r requirements.txt
pip install --index-url https://test.pypi.org/simple/ pele_platform
cd tests
pytest


