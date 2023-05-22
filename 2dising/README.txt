************************************************************************
*** 2dising README file ***
************************************************************************

This is a very simple implementation of the square-lattice Ising model.
We have zero external field (H = 0) and only nearest-neighbor
interactions. All units are reduced units, so J = kB = 1.


COMPILING


You will need a Fortran compiler, as well as the make utility. We recommend
using the GNU compiler, gfortran. No matter your operating system, you will
have to use the command line to run the code. If you are on Ubuntu Linux,
installing all necessary tools will be very easy:

sudo apt install build-essential gfortran

On Windows, it might be easiest to install the Windows 10 Linux sub system with
the latest version of Ubuntu as your Linux flavor.

On macOS, gfortran should be included with xCode (free on the AppStore).

If you run into troubles, feel free to use the discussion board for the
exercises. You can also post how-tos for other systems, such as Microsoft
Visual Studio.

To compile, change into the source code directory and execute

make

This will compile everything. If you make changes to the code, you have
to execute make again to re-compile. To remove the compiled file, use

make clean


USAGE

All the options for the program can be set in the ASCII input file
options.dat. The first line is ignored for comments. Also, the order of
lines is important, and everything after the first item in each line is
ignored as well.

Note that the program does not implement anything to deal with
equilibration. If you want to throw away the initial phase of your
simulation, you will have to either modify the code or read the output,
magnetization.dat, with an external script and do it yourself.

