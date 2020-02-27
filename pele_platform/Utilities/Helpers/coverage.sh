


cd ../../tests

python -m pytest --cov=../ test_others.py  test_defaults.py test_analysis.py -s --cov-report=xml

bash <(curl -s https://codecov.io/bash) -t 8778202a-8b43-4f45-9c6f-82441dc352ec -f coverage.xml
