#!/bin/bash

set -ex

git clone --depth=1 --branch v2.10.0 https://gitlab.dune-project.org/core/dune-common.git || true
git clone --depth=1 --branch v2.10.0 https://gitlab.dune-project.org/core/dune-geometry.git || true
git clone --depth=1 --branch v2.10.0 https://gitlab.dune-project.org/core/dune-grid.git || true
git clone --depth=1 --branch v2.10.0 https://gitlab.dune-project.org/core/dune-localfunctions.git || true
git clone --depth=1 --branch v2.10.0 https://gitlab.dune-project.org/core/dune-istl.git || true
git clone --depth=1 --branch v2.10.0 https://gitlab.dune-project.org/staging/dune-typetree.git || true
git clone --depth=1 --branch v2.10.0 https://gitlab.dune-project.org/staging/dune-uggrid.git || true
git clone --depth=1 --branch v2.10.0 https://gitlab.dune-project.org/staging/dune-functions.git || true
git clone --depth=1 https://gitlab.dune-project.org/extensions/dune-tpmc.git || true

python_minor_version=`python3 -c 'import sys; print(sys.version_info.minor)'`

if [ $? -ne 0 ]; then
	echo 'Executable python3 not available'
	exit 1
fi

if [ $python_minor_version -ge 12 ]; then
	echo 'Python versions starting from 3.12 onwards are not supported'
	echo 'See README for details'
	exit 1
fi

if [ ! -d dune-env ] ; then
	python3 -m virtualenv dune-env
fi

source dune-env/bin/activate
pip install numpy
pip install seaborn
pip install git+https://github.com/tpmc/tpmc.git

cat << EOF > build.opts
CMAKE_FLAGS="-DCMAKE_CXX_COMPILER=g++ \
	     -DCMAKE_C_COMPILER=gcc \
	     -DCMAKE_CXX_FLAGS=\"-O3 -march=native -mtune=native -fno-omit-frame-pointer -g -gdwarf-4 -DHAVE_VARIADIC_TEMPLATES -DNDEBUG -fpermissive -std=c++20\" \
	     -DBUILD_SHARED_LIBS=TRUE \
	     -DTPMC_PREFIX=\"$PWD/dune-env/lib/python3.$python_minor_version/site-packages/tpmc/\""
EOF

./dune-common/bin/dunecontrol --opts=build.opts --builddir=build all
