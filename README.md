# WaveSimulation_Matlab
**Author:** zhaoqingwei  
**Date:** April 24, 2023  
**email:**[zhaoqwei001@163.com](zhaoqwei001@163.com)  


## Main Purpose:
I majored in geophysics with a direction in seismic exploration. Wave Simulation is the key to professional entry, so summarize your own code and share it with you. So the codes is just **simple** **simple** **simple**! The kernel code is just a few lines.
Industrialization code, I recommend C, CUDA, MPI.


## Keyword

* 2-dimension; 3-dimension
* acoustic wave;elastic wave;surface wave
* sponge absorbing boundary condition;Split Perfectly Matched Layer;Convolutional Perfectly Matched Layer
* isotropy,anisotropy media (VTI, TTI)
* Finite Different Modeling 


## Notes

3-dimension program cost memory num*nz*nx*ny*8/1024/1024(MB).Care about memory

while nz=nx=ny=200,num=10; cost memory:610MB
while nz=nx=ny=600,num=50; cost memory:80GB

## Thinks

RToax [https://github.com/Rtoax](https://github.com/Rtoax)  
Guiting-geo [https://github.com/Guiting-geo](https://github.com/Guiting-geo)