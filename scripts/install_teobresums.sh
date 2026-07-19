#!/usr/bin/env bash
set -e

rm -rf build
mkdir build

echo "Downloading TEOBResumS (Dali)..."
git clone -b Dali https://bitbucket.org/teobresums/teobresums.git build/teobresums

# Tested with commit ceb89bc
git checkout ceb89bc

pip install --no-build-isolation ./build/teobresums/Python
pip install --no-build-isolation ./build/teobresums/PyCBC

echo "Successfully installed TEOBResumS"
rm -rf build