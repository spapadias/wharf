#!/bin/bash

# 1. clone yskip repository
cd ../../
git clone https://github.com/yahoojapan/yskip.git
cd yskip

## 2. install yskip
./configure
make
sudo make install

## 3. remove repository
rm -rf ../yskip
