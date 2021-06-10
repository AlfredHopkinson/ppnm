CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt
main: main.o eigenfunc.o
out.txt:main
	./$< > out.txt



Energies.pyxplot.png: energyout.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Energies of the box";\
	plot \
	"energyout.txt" using 1:2 with line title "calculated"\
	,"energyout.txt" using 1:3 with line title "exact"\
	' | tee test.gpi | pyxplot


Modes.pyxplot.png: quantumeigenvalues.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Particle in a box eigen modes";\
	plot \
	"quantumeigenvalues.txt" using 1:2 with line title "Ground State"\
	,"quantumeigenvalues.txt" using 1:3 with line title "First"\
	,"quantumeigenvalues.txt" using 1:4 with line title "Second"\
	' | tee test.gpi | pyxplot




.PHONEY:clean
clean:	
	$(RM) main.o main out.txt

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
