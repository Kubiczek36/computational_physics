************************************************************************
*** lj-gibbs README file ***
************************************************************************


The provided code was originally written by Prof. Martin Neumann and
later slightly modified. It is available for C and Fortran 90.

lj-gibbs performs a Gibbs ensemble Monte Carlo simulation of
Lennard-Jonesium and calculates thermal averages for both volumes. As
the number of particles in each volume is not constant, a trajectory
file is not written, except for the final states when calling agemclj.


COMPILING


You will need a C (or Fortran) compiler, as well as the make utility. We
recommend using the GNU compilers, gcc or gfortran. No matter your
operating system, you will have to use the command line to run the
code. If you are on Ubuntu Linux, installing all necessary tools will be
very easy:

sudo apt install build-essential gcc gfortran

On Windows, the Linux sub system for Windows 10 and a default
installation of Ubuntu are a good option. As of writing this README, the
code does not compile with a default installation of MinGW.

On macOS, gcc should be included with xCode (free on the AppStore).

If you run into troubles, feel free to use the discussion board for the
exercises. You can also post how-tos for other systems, such as
Microsoft Visual Studio.

To compile, change into the source code directory and execute

make

This will compile everything. If you make changes to the code, you have
to execute make again to re-compile. To remove object files, use

make clean

To remove all compiled files:

make realclean

For the Fortran version, only 'make clean' is enough.


USAGE


The simulation code consists of a number of smaller helper programs and
the main Monte Carlo simulation. Usage is the same for both versions,
but binary files are not compatible.

ggemclj
        create a binary input file with a random (gas) configuration

gemclj
        the main simulation code

zgemclj
        "zero out" all accumulated averages after equilibration. Also
        used for changing parameters.

agemclj
        print results in ASCII form, such as accumulated averages

General usage is exactly the same as with the regular Monte Carlo
program, lj-canonical. There are a couple of additional, Gibbs ensemble
specific input parameters, which will be explained here briefly. As with
lj-canonical, reduced LJ units are employed troughout the code. One pass
consists of N = N1 + N2 moves, which can be either single-particle,
volume transfer, or particle transfer moves, chosen at random with a set
probability (see below).

Input parameters (ggemclj):

n1              (initial) number of particles in volume 1
n2              (initial) number of particles in volume 2. Note that in
        GEMC, the total number of particles is fixed.
v1              (initial) volume 1
v2              (initial) volume 2. Note that in GEMC, the total volume
        V = V1 + V2 is fixed.
t               temperature
disp1           maximum displacement for single-particle moves in V1
disp2           maximum displacement for single-particle moves in V2
dv              maximum volume displacement for volume-transfer moves
prob(vmove)     probability to attempt a volume-transfer move
prob(trans)     probability to attempt a particle transfer move
ntskip          Averages are only evaluated every ntskip passes.
ntprint         Current values are printed every ntskip*ntprint passes.
ntjob           total number of passes is ntskip*ntjob

Output parameters (agemclj):

accrpi          Single-particle move acceptance rate in volume i
accrv           volume exchange move acceptance rate
accrt           particle transfer acceptance rate
<Ni>            Average number of particles in volume i
<Vi>            Average volume i
<rho1>          Average density in volume i
<Ui/Ni>         Average potential energy per particle in volume i
<pi>            Average pressure in volume i
<mui>           Average chemical potential in volume i

Some of these values are also printed while gemclj is running. See the
first line for the column labels! Note that the order is slightly
different from the output order in agemclj.

