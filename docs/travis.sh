#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia --color=yes -e 'using Pkg; Pkg.add("Documenter"); Pkg.instantiate(); import ElectromagneticFields; include(joinpath(dirname(pathof(ElectromagneticFields)), "..", "docs", "make.jl"))';
fi
