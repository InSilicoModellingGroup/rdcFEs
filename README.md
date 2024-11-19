# rdcFEs

A simple C++ **r**eaction-**d**iffussion-**c**onvection **F**inite **E**lement **s**olver that utilizes open-source numerical libraries [*PETSc*](https://petsc.org/) and [*libMesh*](https://libmesh.github.io/).
To install these dependencies, the following commands can be executed sequentially, however, the installation directory of *rdcFEs* has to be adapted accordingly.

### Set-up PETSc

First clone the library in some local directory and then checkout to a working version of *PETSc* that works well with *libMesh*.

```
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc && git checkout 746207a3cb30cd796a872832cd0f2ba14a6b454b
```

Subsequently, configure *PETSc* using the command below. Note however, the installation path needs to be provided explicitly as well as the directory of the [*MPI*](https://www.mpich.org/) library intallation.

```
./configure \
  --disable-shared --with-clanguage=C --with-fc=0 \
  --download-cmake --download-f2cblaslapack --download-metis --download-parmetis --download-superlu_dist \
  --with-petsc-arch=macOS__opt \
  --with-mpi-dir=/opt/homebrew/Cellar/open-mpi/4.1.6/ \
  --prefix=/home/USERNAME/petsc-746207a3cb30cd796a872832cd0f2ba14a6b454b
```
It is optional however to provide a custom label for the library set-up (see the '--with-petsc-arch' option above).

Then, follow the commands to run `make` and `make install` that *PETSc* indicates in the terminal after configuration and compilation is completed respectively.

### Set-up libMesh

First clone the library in some local directory and then checkout to a particular version of *libMesh* that *rdcFEs* has been tested to work as intended.

```
git clone https://github.com/libMesh/libmesh.git libmesh
cd libmesh && git checkout d3bda6c7009c6b3241ef2e6c999eff577116dc68 && git submodule update --init --recursive
```

Finally, configure, compile and install *libMesh* using the commands below.

```
./configure \
  --disable-shared --disable-fortran --enable-exodus --enable-metis \
  --with-metis=PETSc --enable-petsc-required \
  PETSC_DIR=/home/USERNAME/petsc-746207a3cb30cd796a872832cd0f2ba14a6b454b \
   --with-cc=/opt/homebrew/Cellar/open-mpi/4.1.6/bin/mpicc  \
  --with-cxx=/opt/homebrew/Cellar/open-mpi/4.1.6/bin/mpicxx \
  --prefix=/home/USERNAME/libmesh-d3bda6c7009c6b3241ef2e6c999eff577116dc68
make -j && make install
```

Do note in the configuration above, however, the explicit provision of the *MPI*-compilers for C and C++ respectively. The user needs to provide the installation path of *libMesh* and indicate where *PETSc* has been installed previously.

### Build rdcFEs

Following installation of the above libraries, *rdcFEs* can be compiled simply by typing `make` in the terminal (in *rdcFEs*' directory of course).
Before that, make sure in the Makefile file the macro `LIBMESH_DIR` is set up correctly.

Then you are ready to go...  🚀
