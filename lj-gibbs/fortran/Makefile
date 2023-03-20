#**********************************************************************#
#
# Makefile for gemclj/f90
#
# 02-Jun-2020 (CL)
#
#**********************************************************************#

F90=gfortran
F90FLAGS=-O3

all: gemclj agemclj ggemclj zgemclj

gemclj: gemclj.f90
	${F90} ${F90FLAGS} -o gemclj gemclj.f90

agemclj: agemclj.f90
	${F90} ${F90FLAGS} -o agemclj agemclj.f90

ggemclj: ggemclj.f90
	${F90} ${F90FLAGS} -o ggemclj ggemclj.f90

zgemclj: zgemclj.f90
	${F90} ${F90FLAGS} -o zgemclj zgemclj.f90
	
clean:
	rm -fv getval_m.mod gemclj_glbm.mod gemclj_subm.mod
	rm -fv gemclj agemclj ggemclj zgemclj

#**********************************************************************#
