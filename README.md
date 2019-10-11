# NUEN-618
Multi-physics simulation
Simple PRKE solver with feedback

Solves the PRKE's via fixed point iteration. Inside the FPI staggered operator splitting or simultaenous solves are used to resolve the coupled nautre of the problem.
Each solve inside FPI requries two non linear solves (one for fuel temperature, one for coolant temperature). Python's scipy.fsolve, gmres, and explicit newtons methods are used.
