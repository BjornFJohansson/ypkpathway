#!/bin/bash


cd $(cd -P -- "$(dirname -- "$0")" && pwd -P)

#export http_proxy='http://proxy.uminho.pt:3128/'
#export https_proxy='https://proxy.uminho.pt:3128/'
#export no_proxy=localhost,127.0.0.0/8,*.local




libreoffice "-env:UserInstallation=file:///tmp/LibO_Conversion" --headless --invisible --convert-to pdf  ./docs/*.odt --outdir ./docs 

libreoffice "-env:UserInstallation=file:///tmp/LibO_Conversion" --headless --invisible --convert-to docx ./docs/*.odt --outdir ./docs

python setup.py register

#python setup.py sdist --formats=gztar,zip
#python setup.py bdist_wheel
#python setup.py bdist_egg
#wine C:/Python27/python.exe setup.py build --plat-name=win32 bdist_wininst
#wine C:/Python27/python.exe setup.py build --plat-name=win-amd64 bdist_wininst


python setup.py sdist --formats=gztar upload
python setup.py sdist --formats=zip upload
python setup.py bdist_egg upload
python setup.py bdist_wheel upload

#wine C://Python27//pythonw.exe setup.py build --plat-name=win32 bdist_wininst upload
#wine C://Python27//pythonw.exe setup.py build --plat-name=win-amd64 bdist_wininst upload


