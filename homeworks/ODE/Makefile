CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

#default: out.txt

#main: main.o func.o 

#out.txt:main
#	./$< > out.txt


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
