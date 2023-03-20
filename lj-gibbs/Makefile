#**********************************************************************#
#
# Makefile for mclj/c
#
# 03-Apr-2010 (MN)
# 19-Mar-2012
#
#**********************************************************************#

.PHONY: clean

CC=gcc
CFLAGS=-O -Wno-unused-result
LFLAGS=-lm

src = $(wildcard *.c)
obj = $(src:.c=.o)
binaries = gemclj ggemclj zgemclj agemclj

all: $(binaries)

$(binaries): %: %.o getval.o
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	rm -f $(obj)

realclean:
	rm -f $(obj) $(binaries)

#**********************************************************************#
