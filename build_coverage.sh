#!/bin/bash -f

if [[ -d BUILD ]]
then
    rm -rf BUILD
fi

export GORILLA_COVERAGE=TRUE

mkdir -p BUILD
cd BUILD
cmake .. -DCMAKE_PREFIX_PATH=$PFUNIT_DIR
make -j

ctest --output-on-failure
cd SRC/CMakeFiles/GORILLA.dir/

gcov-9 *.gcno
lcov --gcov-tool gcov-9 --capture --no-recursion --directory . --output-file covered.info
lcov --gcov-tool gcov-9 --capture --no-recursion -i --directory . --output-file uncovered.info
lcov -a covered.info -a uncovered.info --output-file result.info
genhtml --output-directory ../../../COVERAGE result.info

