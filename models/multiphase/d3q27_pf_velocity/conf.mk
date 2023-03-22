ADJOINT=0
TEST=FALSE
OPT="mapPhi*q27*OutFlow*geometric*BGK*thermo*planarBenchmark*autosym"
# map_phi - maps the hyperbolic tangent function to a linear function for calculating derivatives
# q27 - Q27 lattice structure for phasefield
# ML  - export densities for machine learning
# OutFlow - include extra velocity stencil for outflowing boundaries
# geometric - geometric contact angle implementation, implemented by dmytro sashko
# BGK - single relaxation time operator
# thermo - include energy equation solver for temperature field, influences through
#        - the surface tension
# planarBenchmark - thermocapillary benchmark case
# autosym - symmetry boundary conditions