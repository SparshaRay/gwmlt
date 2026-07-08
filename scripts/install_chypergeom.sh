#!/usr/bin/env bash
set -e

python3 -c '
import juliapkg
juliapkg.require_julia("1.12.6")
juliapkg.add("HypergeometricFunctions", version="0.3.29")
juliapkg.resolve()
juliapkg.status()
'