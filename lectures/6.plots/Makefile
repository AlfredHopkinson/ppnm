CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


gamma.pyxplot.png: gamma.txt Makefile inputgamma.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	"gamma.txt" using 1:2 with line title "tgamma from math.h",\
	"inputgamma.data.txt" using 1:2 with points title "tabulated data";\
	' | pyxplot


erf.pyxplot.png: data.txt Makefile input.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	"$<" using 1:2 with line title "error function"\
	,"$<" using 1:3 with line title "gsl-sf-er"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"input.data.txt" using 1:2 with points title "tabulated data"\
	' | tee test.gpi | pyxplot

test.gpi: erf.pyxplot.png



data.txt : main
	./$< > $@

gamma.txt : main
	./$< > $@

main: main.o erf.o

test.gpi: gamma.pyxplot.png


.PHONEY: clean
clean: 
	$(RM) *.o main*.png data.txt test



test:
	echo $(CFLAGS)
	echo $(LDLIBS)
