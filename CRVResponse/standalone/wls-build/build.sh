#/grid/fermiapp/products/mu2e/artexternals/cmake/v2_8_12_2/Linux64bit+2.6-2.5/bin/cmake -DGeant4_DIR=$GEANT4_FQ_DIR/lib64/Geant4-9.6.3/ ../wls
/grid/fermiapp/products/mu2e/artexternals/cmake/v2_8_4/Linux64bit+2.6-2.5/bin/cmake -DGeant4_DIR=$GEANT4_FQ_DIR/lib64/Geant4-9.6.3/ ../wls

sed 's/\/usr\/local\/lib\/libexpat.so/\/usr\/lib\/libexpat.so/g' CMakeFiles/wls.dir/build.make > CMakeFiles/wls.dir/build.make.tmp
mv CMakeFiles/wls.dir/build.make.tmp CMakeFiles/wls.dir/build.make

echo "CXX_DEFINES += " $1 >> CMakeFiles/wls.dir/flags.make

make
