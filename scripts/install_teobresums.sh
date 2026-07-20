#!/usr/bin/env bash
set -e

rm -rf build
mkdir build

echo "Downloading TEOBResumS (Dali)..."
git clone -b Dali https://bitbucket.org/teobresums/teobresums.git build/teobresums

cd build/teobresums
# Pinned to commit ceb89bc
git -c advice.detachedHead=false checkout ceb89bc

pip install --no-build-isolation ./Python
pip install --no-build-isolation ./PyCBC

echo "Successfully installed TEOBResumS"
rm -rf ../../build