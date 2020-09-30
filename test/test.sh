#!/bin/bash
MAIN_DIR = $pwd
echo 'my dir $MAIN_DIR'
mkdir build && cd build
cmake ..
make all
cd ..
./build/bin/o_tst

