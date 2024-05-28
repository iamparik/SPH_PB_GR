# SPH_PB_GR

##The project is structured as follows:

SPH_PB_GR/
│
├── .gitignore
├── README.md
├── LICENSE
├── CMakeLists.txt (# Optional: if you provide CMake support)
├── depracated
├── docs/ (# Documentation)
├── etc/ (#Directory for configuration files and deployment scripts)
├── examples/ (# Example Cases)
│ └── Case0
│ └── vs_project/ (# Visual Studio project files)
│ ├── Case0.sln
│ └── Case0.vfproj
├── src/ (# Source code)
│ ├── cpp
│ │ ├── mainFile
│ │ ├── headers
│ │ └── functions
│ ├── Fortran
│ │ ├── mainFiles
│ │ ├── modules
│ │ └── subroutines
│ └── python
├── tests/ (# Test files)
│ ├── Periodic
│ │ └── vs_project/ (# Visual Studio project files)
│ │ ├── PeriodicTest.sln
│ │ └── PeriodicTest.vfproj
│ └── FreeSurface
└── ............