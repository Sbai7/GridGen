# GridGen
[![DOI](https://zenodo.org/badge/23898/Sbai7/GridGen.svg)](https://zenodo.org/badge/latestdoi/23898/Sbai7/GridGen)

GridGen is a command-line elliptic grid generation program which is designed to quickly setup logically-rectangular (or IJK ordered) finite volume/element hexahedral meshes for three-dimensional stratigraphic models of porous media subsurface aquifers. It is compatible with GwMove hydrogeological modelling software and could be easily tweaked to adapt to input requirements by other community hydrogeological (or earth-sciences) computer modelling codes. 

The algorithm is based on numerical solution of an elliptic PDE namely the Laplace 2D equation in the plane [1,2,3]. Next, to build the full 3D grid, the planar grid is mapped into intermediate slices separating hydrogeological layers with different lithology. Intermediate slices in each hydrogeological unit may be specified to get required accuracy along the vertical direction. The method is proven for its ability to generate very smooth and very high quality meshes for demanding finite element solvers. 

The power of GridGen lies in simplicity of usage and lightweight portable code which could be integrated easily into the core of legacy Fortran programs for automated boundary-fitted grid generation. However, it's not a universal automatic mesh generator since you cannot mesh general domains especially those with internal geometric constraints. Even with these limitations the program has proven suitable for many groundwater modelling projects and has been used extensively by the authors and their colleagues. 

GridGen is brought to you to be useful by helping you not to reinvent the wheel. We'll be very pleased to use the following citation when reporting results with this program:

Sbai, M.A. (2016). GridGen v1.0. Zenodo. 10.5281/zenodo.60439 

## Building GridGen from source 

### Supported Platforms / Compilers 
The program has been tested on the following platforms:
- Microsoft Windows XP/7/10 
- Ubuntu Linux 
- Cygwin (Linux emulation layer for Windows platforms)

The following compilers are equally supported: 
- Intel Visual Fortran 11.1 or above 
- GNU Fortran 4.2 or above 

However only the following Platform/compiler combination has been tested to date:
- Windows + Intel Visual Fortran 
- Ubuntu Linux + GNU Fortran
- Cygwin + GNU Fortran 

Follow the instructions below to compile GridGen in your platform.

### Windows Platforms 
You must have installed Visual Studio 2005 or above with a compatible Intel Visual Fortran compiler. You must select the installation of the Visual Studio integration package during the last installation to build VS Fortran projects right away from the graphical user interface. 
- Double click on the Visual Studio 2005 solution file '*GridGen.sln*' to launch it in Visual Studio 2005. For later versions, accept the project upgrade process and save your newly converted solution.
- Select '*Release*' as a solution configuration and '*Win32*' as a solution platform from the two dropdown menus in the main toolbar. 
- Fire 'F7' key or select '*Build > Build Solution*' menu command, then wait for a few minutes until compilation of all the project files. 
- The target executable is located in '*bin*' subdirectory of the distribution root folder. 

That's all!

### Linux / Cygwin Platforms 
Go to the root folder first, then copy & paste the following sequence of commands to your console window or a new shell file created by your own:

```
mkdir build 
cd build
gfortran -c ../src/string.f90 
gfortran -c ../src/globals.f90
gfortran -c ../src/bndFitted.f90
gfortran -c ../src/command_line_syntax.f90
gfortran -c ../src/main.f90
gfortran -c ../src/parse_command_line.f90
gfortran -c ../src/writeMesh.f90
gfortran -c ../src/xybord.f90
gfortran -c ../src/xyfixed.f90
gfortran -c ../src/xyfree.f90
gfortran -c ../src/zfixed.f90
gfortran -c ../src/zfree.f90
gfortran *.o -o gridgen
rm *.o *.mod 
```

Finally, install the binary file into any folder in your system path (e.g. ``` ~/bin ```):

```
mv gridgen ~/bin 
```

## Binary Releases
We offer binary releases for windows platforms (>7) and they may be found in the '*releases*' page of this github repository. 

## Runing GridGen 

### Input file 

### Output file(s) 

## Examples 

### 3D grid of a soil sample

Here is an example of grid construction fitting the boundaries of a soil sample. First, we start by 2D gridding in the plane where 20 points of the circular circumference have been specified among which only four corner points have been explicitly fixed. Finally, the 2D mesh is simply projected in the vertical direction with equidistant spacing. Here it is:

![Alt text](pictures/example1.jpg?raw=true "")

This is not common in common hydrogeological applications but it is valuable in other contexts such as for two-phase / multiphase flow simulations to replicate laboratory experiments with other dedicated engines. We'll post a computational example later on. 

### 3D radial grid arround a wellbore 

Another, commonly encoutered situation in hydrogeology is related to well testing of hydrogeological characteristics of single and/or multiple aquifer systems. Analytical solutions for special cases are standard practice. To go even further by alleviating all the limitations of available analytical solutions a numerical groundwater flow model is more appropriate. In this context GridGen is also useful since it has the ability to generate high qualitiy grids arround the wellbore. The area in the the internal well radius and a given far flow boundary at which the drawdown is negligeable is gridded. The radial spacing could follow a linear, logarithmic or geometric distribution as shown in the following figure:


### 3D grid of a realistic multilayer aquifer 

This example demonstrates typical use of the program in the context of a realistic application designed to model flow and solute transport in a multilayer aquifer. The upper unit belongs to an uncofined aquifer whose mean depth is 11.7m while the lower confined aquifer mean depth is 5.7m. The two aquifers are separated by an aquitard with uniform thickness of 4.5m. A landfill contaminant area (not shown here) is located at the upper slice. The 'snapy nodes' features of the code was used to snap the closest nodes to this inner boundary and locate all nodes inside this polygon. The figure showing the aquifer geometry is just diplayed below:

![Alt text](pictures/example4_geom.jpg?raw=true "")

Here is a closeup view of the full 3D mesh which was designed for the groundwater simulation. It is composed from 73 nodes (or 72 cells), 89 nodes (or 88 cells) and 12 slices (11 layers) aligned along the X, Y and Z directions respectively. Thus, forming a mesh whose total number of nodes is 77964 and total number of hexahedral cells is 69696. 

![Alt text](pictures/example4_mesh.jpg?raw=true "")

Here is an example of the steady-state groundwater flow simulation with GwMove software computed on the previously constructued grid. Details of the simulation will be documented in GwMove wiki pages. 

![Alt text](pictures/example4_heads.jpg?raw=true "")

## Acknowledgments 
The initial F77 program (Geo_Grid) was developed by Prof. Abdelkader Larabi (now @ Ecole Mohammedia d'Ingénieurs, University M. V, Rabat, Morocco) during his Ph.D. Thesis between 1990 and 1994 at the laboratory of Hydrology, Free University, Brussels. It was improved by Dr. M. Adil Sbaï (now @ French Geological Survey) during his Ph.D. thesis days between 1995 and 1999 in the same department. Since then the program was used extensively in many groundwater modelling projects. 

## References 

[1] Thompson, J.F.; Thames, F.C.; Mastin, C.W. (1974). "Automatic Numerical Generation of Body-fitted Curvilinear Coordinate System for Field Containing any Number of Arbitrary Two-Dimensional Bodies". J. Comput. Phys. 15: 299–319. doi:10.1016/0021-9991(74)90114-4.

[2] Thompson, Joe F., Warsi, Z. U. A. and Mastin, C. W. (1982). "Boundary-Fitted Coordinate Systems for Numerical Solution of Partial Differential Equations -- A Review", J. Comput. Phys., 47(1). 

[3] Thompson, Joe F. (Ed.) Numerical Grid Generation, North-Holland 1982. (Also published as Vol. 10 11 of Applied Mathematics and Computation, 1982).
