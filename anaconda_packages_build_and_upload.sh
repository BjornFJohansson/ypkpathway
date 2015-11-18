#!/bin/bash

echo "Building anaconda packages"

conda build .

OUTPUT="$(conda build . --output)"
echo "${OUTPUT}"

conda convert $OUTPUT -p linux-32 --output-dir ~/anaconda/conda-bld/
conda convert $OUTPUT -p win-32   --output-dir ~/anaconda/conda-bld/
conda convert $OUTPUT -p win-64   --output-dir ~/anaconda/conda-bld/
conda convert $OUTPUT -p osx-64   --output-dir ~/anaconda/conda-bld/


exit 42


anaconda upload $OUTPUT

name=$(basename $OUTPUT)

anaconda upload ~/anaconda/conda-bld/linux-32/$name
anaconda upload ~/anaconda/conda-bld/win-32/$name
anaconda upload ~/anaconda/conda-bld/win-64/$name
anaconda upload ~/anaconda/conda-bld/osx-64/$name

# $SHELL



# http://docs.binstar.org/conda.html
# If you have previously generated TOKEN (check Token Generation) then you may run upload process also in this way:
# $ binstar -t ${TOKEN} upload /home/USERNAME/anaconda/conda-bld/linux64/conda_gc_test-1.2.1-py27_3.tar.bz2
# http://docs.binstar.org/conda.html
