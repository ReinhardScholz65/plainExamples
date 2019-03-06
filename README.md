## Installation

* Clone and Install the library `git clone https://github.com/LeMonADE-project/LeMonADE.git`
* Install cmake (minimum version 2.8)
* Go to one of the project directories and perform the standard build:
 
````sh
    # prepare the build directory
    mkdir build
    cd build
    # run cmake with the lemonade path
    cmake -DLEMONADE_DIR=</path/to/LeMonADE-library/>  ..
    make
````

## Examples

* in analyzers and updaters examples files can be found
* the are used in [the example project](projects/exampleProject/)

## Tests

* in [the test directory](test/) you find a setup using the [catch2](https://github.com/catchorg/Catch2) framework 
* run the tests similar to the projects
````sh
    # prepare the build directory
    mkdir build
    cd build
    # run cmake with the lemonade path
    cmake -DLEMONADE_DIR=</path/to/LeMonADE-library/>  ..
    make
    ./testsYves
````
* add tests by creating a **test_yourtestname.cpp** file and add this to the [CMakeLists.txt](test/CMakeLists.txt):

```cmake
    add_executable(testsLeMonADE test_main.cpp .... )
```

* testing might be annoying but helps you to find nasty bugs
 

