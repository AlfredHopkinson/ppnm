CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

#default: out.txt

#main: main.o 

#out.txt:main
#	./$< > out.txt


#.PHONEY:clean
#clean:	
#	$(RM) main.o main out.txt simple_out.txt 


Simple.pyxplot.png: out.txt simple_out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Example of the subspline on a simple generated points";\
	plot \
	"simple_out.txt" using 1:2 with line title "Akima sub-spline"\
	,"simple_points.txt" using 1:2 with points title "given points"\
	' | tee test.gpi | pyxplot




test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
