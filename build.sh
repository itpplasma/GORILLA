#!/bin/bash -f

if [[ -d BUILD ]]
then
    rm -rf BUILD
fi

export GORILLA_COVERAGE=FALSE

mkdir -p BUILD
cd BUILD
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j

