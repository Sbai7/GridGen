# GridGen
[![DOI](https://zenodo.org/badge/23898/Sbai7/GridGen.svg)](https://zenodo.org/badge/latestdoi/23898/Sbai7/GridGen)

GridGen is a command-line grid generation program which is designed to quickly setup logically-rectangular (or IJK ordered) finite volume/element meshes for three-dimensional stratigraphic models of pourous media subsurface aquifers. It is compatible with GwMove hydrogeological modelling software and could be easily tweaked to adapt to input requirements by other community hydrogeological (or earth-siences) computer modelling codes. 

Please use the following citation when reporting results with this software:

Sbai, M.A. (2016). GridGen v1.0. Zenodo. 10.5281/zenodo.60439 

## Building GridGen from source 

### Supported Platforms / Compilers 
The program has been tested on the following platforms:
- Microsoft Windows XP/7/10 
- Ubuntu Linux 

The following compilers are equally supported: 
- Intel Visual Fortran 11.1 or above 
- GNU Fortran 4.2 or above 

However only the following Platform/compiler combination has been tested to date:
- Windows + Intel Visual Fortran 
- Ubuntu Linux + GNU Fortran

Follow the instructions below to compile GridGen in your platform.

### Windows Platforms 
You must have installed Visual Studio 2005 or above with a compatible Intel Visual Fortran compiler. You must select the installation of the Visual Studio integration package during the last installation to build VS Fortran projects right away from the graphical user interface. 
- Double click on the Visual Studio 2005 solution file '*GridGen.sln*' to launch it in Visual Studio 2005. For other version, accept the project upgrade process and save your new solution.
- Select '*Release*' as a solution configuration and '*Win32*' as a solution platform from the two dropdown menus in the main toolbar. 
- Fire 'F7' key or select '*Build > Build Solution*' menu command, then wait for a few minutes until compilation of all the project files. 
- The target executable is located in '*bin*' subdirectory of the git distribution root folder. 

That's all !

### Linux Platforms 
- 
- 
- 
- 

## Binary Releases
We offer binary releases for windows platforms (>7) and they may be found in the '*releases*' page of this github repository. 

## Runing GridGen 

### Input file 

### Output file(s) 

## Examples 