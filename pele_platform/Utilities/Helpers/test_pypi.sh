#Clean after build
pip uninstall pele_platform --yes
pip install --index-url https://test.pypi.org/simple/ pele_platform
cd tests
pytest


