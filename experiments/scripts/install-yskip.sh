#!/bin/bash

# clone yskip repository
cd ../
git clone https://github.com/yahoojapan/yskip.git
cd yskip

# install yskip
./configure
make
sudo make install

# remove repository
rm -rf ../yskip