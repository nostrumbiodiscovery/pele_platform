rm -rf AdaptivePELE
git clone https://github.com/AdaptivePELE/AdaptivePELE.git
cd AdaptivePELE
conda build -c nostrumbiodiscovery -c conda-forge --python 3.8 conda_recipe
cd ..

rm -rf PlopRotTemp
git clone https://github.com/NostrumBioDiscovery/PlopRotTemp.git
cd PlopRotTemp
conda build  conda_recipe/ --python 3.8
cd ..

rm -rf ProDy
git clone https://github.com/danielSoler93/ProDy.git
cd ProDy
conda build conda_recipe/ --python 3.8
cd ..

rm -rf PPP
git clone https://github.com/NostrumBioDiscovery/PPP.git
cd PPP
conda build -c nostrumbiodiscovery -c conda-forge conda_recipe/ --python 3.8
cd ..


rm -rf pyfpdf
git clone https://github.com/danielSoler93/pyfpdf.git
cd pyfpdf
conda build conda_recipe/ --python 3.8
cd ..

conda build -c nostrumbiodiscovery -c conda-forge conda_recipe/ --python 3.8 

