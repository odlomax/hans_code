# hans_code
Some code for my buddy, Hans.

In order to build the Fortran part of this code, you'll need a Fortran compiler.

Before you run the Python code (interpolate_test.py), you'll need to execute the command:

python -m numpy.f2py -c f_module.pyf interpolate.f90
