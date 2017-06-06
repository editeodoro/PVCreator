# PVCreator
Create slices of a datacube

Just set the cfitsio and wcslib path in the Makefile. Then:

> make 

> make install (if you want to copy the exe somewhere)

Enjoy.


PVCREATOR usage: 

   1) pvcreator <inp_cube> <out_pv> <x0> <y0> <angle> 
   2) pvcreator <inp_cube> <out_pv> <x1> <y1> <x2> <y2> 

In 1) the slice is defined by the point (x0,y0) and angle (N->W).
In 2) the slice is defined by two points (x1,y1) and (x2,y2)



