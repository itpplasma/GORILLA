if [[ -d BUILD ]]
then
    rm -rf BUILD
fi

export GORILLA_COVERAGE=FALSE

mkdir -p BUILD
cd BUILD
cmake ..
make -j

