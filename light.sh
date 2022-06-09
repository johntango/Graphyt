#! /bin/bash
cmake .
make -j12
python -m unittest discover ./test/ -v
