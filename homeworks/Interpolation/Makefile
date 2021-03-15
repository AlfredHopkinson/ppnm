CC = gcc
LDLIBS = -lm
CFLAGS = -Wall -O1 -std=gnu11


#linear.out.txt: linear
#	./$< > $@


#points.pyxplot.png: plots.out.txt Makefile
#	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Random Points";\
	plot \
	"plots.out.txt" with points using 1:2\
	' | tee test.gpi | pyxplot

#interpolation.pyxplot.png: plots.out.txt linear.out.txt Makefile
#	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Interpolation";\
	plot \
	"plots.out.txt" using 1:2 with points\
	,"linear.out.txt" with line title "linear interpolation"\
	' | tee test.gpi | pyxplot


integration.pyxplot.png: linearint.out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Linear Integration";\
	plot \
	"linearint.out.txt" with line using 1:2\
	' | tee test.gpi | pyxplot

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
