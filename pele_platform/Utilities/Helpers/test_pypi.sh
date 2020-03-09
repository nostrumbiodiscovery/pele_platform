#Clean after build
conda create -n test_pele_plat python=3.7 --yes
activate test_pele_plat
pip uninstall -r requirements.txt --yes
pip uninstall pele_platform --yes
pip install numpy cython
pip install .
cd tests
pytest
conda remove --name myenv 

