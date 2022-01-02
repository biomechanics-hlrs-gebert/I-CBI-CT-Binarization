# I-CTBI-CT_Binarization
![Architecture](https://img.shields.io/badge/Architecture-x86-green)
![OS](https://img.shields.io/badge/Linux-64Bit-green)
![version](https://img.shields.io/badge/version-1.0.0-green)
![Contributors](https://img.shields.io/badge/HLRS-NUM-blue)

This program reads 3-Dimensional scalar fields to binarize and eventually invert them.

It's tested with up to 160 Processors (4 nodes) on Vulcan. Turnaround time of roughly 200 Seconds while reading/computing/writing to storage on 18.5E09 Voxels of kind INTEGER2.

The program currently support integer kind=2 and integer kind=4 (ik2/ik4) scalar fields.
## Meta Template
Located in: 
```
./datasets/I-CTBI.meta.template
```
For use with previously used data sets:
```
cat ./datasets/I-CTBI.meta.template >> Your_Meta_File.meta
```

## [semantic versioning](https://semver.org):

Given a version number MAJOR.MINOR.PATCH, increment the:

* MAJOR version when you major Features (i.e. new way of image processing),
* MINOR version when you extend functionality (i.e. new kernels), and
* PATCH version when you make bug fixes.

## Requirements
* x86 64bit Hardware
* Linux x86 64Bit Installation with a BASh
* GNU Compiler Collection (GCC) including gcc/gfortran
* An installation of Open MPI - run the script in the project's root directory.
### Message Passing Interface 
Parallelization of the program is done with an API called MPI (Message Passing Interface).

Required: MPI - compiled with integer 4 and mpi_f08.

  1. [Open-mpi 4.1.0](https://www.open-mpi.org/software/ompi/v4.1/) on local systems. Other versions are not tested.
  2. [HPE-MPT on HLRS Hawk](https://kb.hlrs.de/platforms/index.php/MPI(Hawk))

The program may be ported to other architectures. Maybe not :-)

## Build
It's tested and therefore recommended to build and run the program as follows.
### Set up the Environment
```vim ./auxiliaries/system_environments/<system>.sh```
```source ./environment.sh <system>``` 

* Set an architecture/a system
  * Give the absolute base path of your mpi-installation
  * Alternatively give the proper module names of your compute cluster

### Run make:
Build the program:    ```make```
Create documentation: ```make docs```

### Uninstall:
```make clean && rm -r <your program directory>```
