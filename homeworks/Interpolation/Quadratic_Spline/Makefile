CC = gcc
LDLIBS = -lm
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
 

quad.out.txt: quadratic
	./$< > $@



interpolation.pyxplot.png: quad.out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Quadratic Interpolation";\
	plot \
	"plots.out.txt" using 1:2 with points\
	,"linear.out.txt" with line title "linear interpolation"\
	,"gslinfo.txt" using 1:2 with line title "linear gsl interpolation"\
	' | tee test.gpi | pyxplot
