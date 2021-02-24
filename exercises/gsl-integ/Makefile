CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

#error.pyxplot.png: out.error.txt Makefile
#	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function intergration";\
	plot \
	"out.error.txt" using 1:2 with line title "Error function";\
	' | pyxplot

bessel.pyxplot.png: out.bessel.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Bessel function intergration";\
	plot \
	"out.bessel.txt" using 1:2 with line title "Bessel function";\
	' | pyxplot



.PHONEY: default
	default: out.main.txt out.params.txt out.mygamma.png
out.main.txt:main
	./$< > $@
out.params.txt:params
	./$< > $@
out.error.txt : error
	./$< > $@
out.bessel.txt :bessel
	./$< > $@


.PHONEY : clean
clean:
	$(RM) out* main params mygamma

test :
	echo $(CFLAGS)
