#!/bin/sh

Geant4_DIR=~/geant4.10.6/lib/Geant4-10.6.1/
CMAKE_MODULE_PATH=$Geant4_DIR/Modules
CMAKE_PREFIX_PATH=/usr/lib/x86_64-linux-gnu/qt4

echo $1

cmake $1 -DGeant4_DIR=$Geant4_DIR -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH -DCMAKE_MODULE_PATH=$CMAKE_MODULE_PATH