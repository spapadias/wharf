#!/bin/bash

# clone yskip repo
cd ../
git clone https://github.com/yahoojapan/yskip.git
cd yskip

# install yskip
./configure
make
sudo make install

# remove repo
rm -rf ../yskip