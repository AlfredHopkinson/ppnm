FLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


gamma.pyxplot.png: gamma.txt Makefile inputgamma.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma Function";\
	plot \
	"gamma.txt" using 1:2 with line title "tgamma"\
	,"gamma.txt" using 1:4 with line title "log(tgamma)"\
	,"inputgamma.data.txt" using 1:2 with points title "tabulated data"\
	' | tee test.gpi | pyxplot

gamma.txt: main
	./$< > $@


.PHONEY: clean
clean: 
	$(RM) *.o main*.png data.txt test



test:
	echo $(CFLAGS)
	echo $(LDLIBS)
