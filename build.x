if [[ -d BUILD ]]
then
    rm -rf BUILD
fi

mkdir -p BUILD
cd BUILD
cmake .. -DCMAKE_PREFIX_PATH=$PFUNIT_DIR
make -j

cd BUILD
ctest --verbose
cd SRC/CMakeFiles/GORILLA.dir/

/usr/bin/gcov-9 *.gcno
lcov --gcov-tool /usr/bin/gcov-9 --capture --directory . --output-file covered.info
lcov --gcov-tool /usr/bin/gcov-9 --capture -i --directory . --output-file uncovered.info
lcov -a covered.info -a uncovered.info --output-file result.info
genhtml --output-directory html result.info


