CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt

main: main.o leastsf.o

out.txt:main
	./$< > out.txt

#LSF.pyxplot.png: LSF.out.txt plot.out.txt
#	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "Time in Days";\
	set ylabel "Log(Activity)";\
	set title "Least Square Fit";\
	plot \
	"LSF.out.txt" using 1:2 with line title "best fit"\
	,"plot.out.txt" using 1:2:3 with errorbars title "data"\
	' | tee test.gpi | pyxplot









.PHONEY:clean
clean:
	$(RM) main.o main out.txt

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
