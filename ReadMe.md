## Build rdcFE

Before building the main code, you need to download, build and install the PETSc and libMesh library

Install these libraries in a separate folder from the code source (set the directory path in the commands found below)

>> Make sure that appropriate verions of gcc, cmake and openmpi are available. In the case of the Turing cluster, the module load commands are found below.

### Compiler versions

These are working on Turing cluster but may not be the only ones

gcc (GCC) 8.3.0
mpirun (Open MPI) 3.1.4
cmake version 3.15.3

### Turing cluster packages

module use /share/apps/eb/modules/all
module use /home/eioann18/.local/easybuild/modules/all
module load EasyBuild

module load foss/2019b
module load CMake/3.15.3-GCCcore-8.3.0
module load mc/4.8.13-GCC-4.9.2

### Configuration for PetSc

```
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc && git checkout 746207a3cb30cd796a872832cd0f2ba14a6b454b
./configure \
  --disable-shared --with-clanguage=C --with-fc=0 --download-cmake --download-f2cblaslapack --download-metis --download-parmetis --download-superlu_dist \
  --with-petsc-arch=centos__opt \
  --with-mpi-dir=/share/apps/eb/software/OpenMPI/3.1.4-GCC-8.3.0/ \
  --prefix=/home/eioann18/repository/libs/petsc_installation
make PETSC_DIR=/home/eioann18/repository/libs/petsc PETSC_ARCH=centos__opt all
make PETSC_DIR=/home/eioann18/repository/libs/petsc PETSC_ARCH=centos__opt install
make PETSC_DIR=/home/eioann18/repository/libs/petsc PETSC_ARCH=centos__opt check
```

### Configuration for libMesh

```
git clone https://github.com/libMesh/libmesh.git libmesh
cd libmesh && git checkout d3bda6c7009c6b3241ef2e6c999eff577116dc68 && git submodule update --init --recursive
./configure \
  --disable-shared --disable-fortran --enable-exodus --enable-metis --with-metis=PETSc --enable-petsc-required \
  PETSC_DIR=/home/eioann18/repository/libs/petsc_installation  \
  --with-cxx=/share/apps/eb/software/OpenMPI/3.1.4-GCC-8.3.0/bin/mpicxx \
   --with-cc=/share/apps/eb/software/OpenMPI/3.1.4-GCC-8.3.0/bin/mpicc  \
  --prefix=/home/eioann18/repository/libs/libmesh_installation 
make && make install
```

### Make rcdFE

Clone for repository
`git clone https://github.com/vasvav/rdcFEs.git`

Set libmesh installation path in Makefile and run
`make all`