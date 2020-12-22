# PIC/FLIP Simulation

Implemented by Yongzhen Huang and Gongyi Shi.

Based on the paper [*Animating Sand as a Fluid*](https://www.cs.ubc.ca/~rbridson/docs/zhu-siggraph05-sandfluid.pdf).

A video presenting our work can be found [here](https://www.dropbox.com/s/6ewmlhst6lxcsgc/FLUID_SIM_PRESENTATION.mp4?dl=0).


## How to build
This project requires [libigl](https://github.com/libigl/libigl) to work.

Clone the repository and download Libigl
```
// go into project folder
git clone --recursive https://github.com/libigl/libigl.git
```

CMake is use to build the project
```
// go into project folder
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
```

## Running the program
Simply execute the program 
```
./Fluid_Sim [path to triangle mesh]
```

Optional parameter of path to a triangle mesh shows the object falling into the fluid.

Program options are shown immediately afterwards.