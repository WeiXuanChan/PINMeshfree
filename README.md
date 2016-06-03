# Point Interpolated Newton-Rapshon Meshfree Method

Keywords: Meshfree, Finite Element Analysis, Computational fluid dynamics

###Requirement

Numpy

AutoD - https://github.com/WeiXuanChan/autoD

#####Basic Requirement (standard module)
-pinm module: sys, matplotlib, time, scipy

-discretizer module: stl, mpl_toolkits, matplotlib

###Description

General numerical solver using automatic differentiation and Point Interpolation Method with Meshfree Method.

Enables user defined governing equation (non-linear accepted), material properties and basis.

###Work flow
Import -> Mesh -> Filter and Link Nodes -> setup Solver (variables, equations, material) -> see results (detailed postprocessing currently in progress)
###Function
######pinm

