CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: plots.out.txt
plots.out.txt: linear
	./%< > $@

interpolation.pyxplot.png: Makefile plots.out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Interpolation";\
	plot \
	"plots.out.txt" using 1:2\
	' | tee test.gpi | pyxplot

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
