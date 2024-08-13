cd /app

source pythonenv3.8/bin/activate
python setup.py bdist_wheel
deactivate

source pythonenv3.9/bin/activate
python setup.py bdist_wheel
deactivate

source pythonenv3.10/bin/activate
python setup.py bdist_wheel
deactivate

source pythonenv3.11/bin/activate
python setup.py bdist_wheel
deactivate

source pythonenv3.12/bin/activate
python setup.py bdist_wheel
deactivate


# auditwheel repair --plat manylinux2014_x86_64 dist/*