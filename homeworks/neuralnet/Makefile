CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

#default: out.txt num.txt

#main: main.o ann.o mini.o

#out.txt:main
#	./$< > out.txt


#.PHONEY:clean
#clean:	
#	$(RM) main.o main out.txt num.txt mini.o mini ann.o ann


NN.pyxplot.png: out.txt num.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "SIR Model";\
	plot \
	"out.txt" using 1:2 with linespoints title "generated points"\
	,"out.txt" using 1:3 with linespoints title "trained"\
	' | tee test.gpi | pyxplot




test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
