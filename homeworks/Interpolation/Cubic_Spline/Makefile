CC = gcc
LDLIBS = -lm
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
 

#cubic.out.txt: cubic
#	./$< > $@


cubic_interpolation.pyxplot.png: plots.out.txt cub.out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Cubic Interpolation";\
	plot \
	"plots.out.txt" using 1:2 with points\
	,"cub.out.txt" using 1:2 with line title "cubic interpolation"\
	,"cub.out.txt" using 1:3 with line title "cubic interpolation differentiation"\
	,"cub.out.txt" using 1:4 with line title "cubic interpolation integration"\
	' | tee test.gpi | pyxplot
