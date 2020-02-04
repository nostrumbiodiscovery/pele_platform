#Clean after build
pip uninstall -r requirements.txt --yes
pip uninstall pele_platform --yes
pip install numpy cython
pip install -r requirements.txt --no-cache-dir
pip install --index-url https://test.pypi.org/simple/ pele_platform --no-cache-dir
cd tests
pytest


