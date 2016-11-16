cd argodsm
cd tests
wget https://github.com/google/googletest/releases/tag/release-1.7.0.zip
unzip release-1.7.0.zip
mv googletest-release-1.7.0 gtest-1.7.0
cd ..
mkdir build
cd build
cmake -DARGO_BACKEND_MPI=ON        \
      -DARGO_BACKEND_SINGLENODE=ON \
      -DARGO_TESTS=ON              \
      -DBUILD_DOCUMENTATION=ON     \
      -DCMAKE_CXX_COMPILER=mpic++  \
      -DCMAKE_C_COMPILER=mpicc     \
      ../
make
make test
