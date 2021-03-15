CC = gcc
LDLIBS = -lm
CFLAGS = -Wall -O1 -std=gnu11


#plots.out.txt: linear
#	./$< > $@


interpolation.pyxplot.png: plots.out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Interpolation";\
	plot \
	"plots.out.txt" using 1:2\
	' | tee test.gpi | pyxplot




test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
