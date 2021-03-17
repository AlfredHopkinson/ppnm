CC = gcc
LDLIBS = -lm
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
 

#quad.out.txt: quadratic
#	./$< > $@



interpolation.pyxplot.png: plots.out.txt quad.out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Quadratic Interpolation";\
	plot \
	"plots.out.txt" using 1:2 with points\
	,"quad.out.txt" using 1:2 with line title "quadratic interpolation"\
	,"quad.out.txt" using 1:3 with line title "quadratic interpolation differentiation"\
	,"quad.out.txt" using 1:4 with line title "quadratic interpolation integration"\
	' | tee test.gpi | pyxplot
