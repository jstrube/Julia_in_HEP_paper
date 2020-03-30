# Julia code for the paper "Performance of Julia for High Energy Physics Analyses"
https://arxiv.org/abs/2003.11952

## Installation
The directory is self-contained and lists all dependencies necessary to run the code.
* Download the latest Julia binary for your platform from https://julialang.org/downloads/
* checkout the code from this repository
* `cd /into/the/directory/with/this/README.md`
* `julia`
* `] activate .`
* `instantiate`
* Exit the julia REPL by pressing Ctrl+D

That should be it. The steps above download all dependencies necessary to run the code, so they only need to be executed once.

The code is then executed using `julia <nameofscript.jl> <lciofile1> [<lciofile2> ...]` 
Note that the first time the code is executed, julia pre-compiles the dependencies. This is done transparently, but it means that the first run will be slower than subsequent runs. 
