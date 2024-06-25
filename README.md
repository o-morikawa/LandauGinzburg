# LandauGinzburg

See [Doxygen HTML](https://o-morikawa.github.io/LandauGinzburg/src/html/index.html).

## Numerical study of the N=2 Landau--Ginzburg models

author Okuto Morikawa

date   Created on May 8 2018

## User's guide
Code set for the generation of the configuration of $N(p)$, and the computation of $A(p)$ and sign determinant

### Environment
Ubuntu 16.04: landau [Xeon E5 2660V4 2.0GHz 28C/56T], lifshitz2\&3 [Xeon Silver 4114 2.20GHz 10C/20T], and supercomputer ITO.

### Command
```g++ -std=c++11 (-O2) -fopenmp```

### Files
- field_class.hpp; field_class.cpp
- field_potential.hpp; field_potential.cpp
- field_nicolai.hpp; field_nicolai.cpp, field_spt.cpp
- field_conf.cpp

### Eigen
These codes are based on [Eigen v3.3.2](http://eigen.tuxfamily.org/).
Please download the Eigen code set as "Eigen" directory in src/src_imp.

### Julia
Computation of some correlation functions by [julia v0.5.2](https://julialang.org/)
