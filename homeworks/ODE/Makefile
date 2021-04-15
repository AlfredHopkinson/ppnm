CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt sir.out.txt

main: main.o func.o  

out.txt:main
	./$< > out.txt

sir: sir.o func.o

sir.out.txt:sir
	./$< > sir.out.txt

#sir_T_small.pyxplot.png: sir.out.txt 
#	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "SIR Model with Tc of 1";\
	plot \
	"sir.out.txt" using 1:2 with linespoints title "Susceptible"\
	,"sir.out.txt" using 1:3 with linespoints title "Infected"\
	,"sir.out.txt" using 1:4 with linespoints title "Removed"\
	' | tee test.gpi | pyxplot




sir.pyxplot.png: sir.out.txt 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "SIR Model";\
	plot \
	"sir.out.txt" using 1:2 with linespoints title "Susceptible"\
	,"sir.out.txt" using 1:3 with linespoints title "Infected"\
	,"sir.out.txt" using 1:4 with linespoints title "Removed"\
	' | tee test.gpi | pyxplot


example.pyxplot.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Example";\
	plot \
	"out.txt" using 1:2 with line\
	' | tee test.gpi | pyxplot



.PHONEY:clean
clean:	
	$(RM) main.o main out.txt

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
