CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

report.pdf: report.tex 
	pdflatex $<
	pdflatex $<
	bibtex $<
	pdflatex $<

cos.pdf: cos.gpi cos_out.txt
	pyxplot $<

intder.pdf: intder.gpi
	pyxplot $<



default: out.txt
	
main: main.o 

out.txt:main
	./$< > out.txt


.PHONEY:clean
clean:	
	$(RM) main.o main out.txt simple_points.txt simple_out.txt cospoint_out.txt cos_out.txt 



Simple.pyxplot.png: out.txt simple_out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Example of the subspline on a simple generated points";\
	plot \
	"simple_out.txt" using 1:2 with line title "Akima sub-spline"\
	,"simple_points.txt" using 1:2 with points title "given points"\
	,"simple_out.txt" using 1:3 with line title "cubic sub-spline"\
	' | tee test.gpi | pyxplot


cos.pyxplot.png: out.txt cos_out.txt cospoint_out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Example of the subspline for cos";\
	plot \
	"cospoint_out.txt" using 1:2 with points title "generated points"\
	,"cospoint_out.txt" using 1:3 with points title "derivitive"\
	,"cos_out.txt" using 1:2 with line title "spline"\
	,"cos_out.txt" using 1:3 with line title "akima"\
	' | tee test.gpi | pyxplot


deriv.pyxplot.png: out.txt intder.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Example of the subspline for cos";\
	plot \
	"intder.txt" using 1:2 with line title "Derivitive"\
	' | tee test.gpi | pyxplot

integ.pyxplot.png: out.txt intder.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Example of the integral of the subspline for cos";\
	plot \
	"intder.txt" using 1:3 with line title "Integral"\
	' | tee test.gpi | pyxplot

mix.pyxplot.png: intder.txt cos_out.txt cospoint_out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "Example of the subspline and its derivitive and integral";\
	plot \
	"cospoint_out.txt" using 1:2 with points title "generated points"\
	,"cos_out.txt" using 1:2 with line lw 2 title "spline"\
	,"intder.txt" using 1:2 with line title "derivitive"\
	,"intder.txt" using 1:3 with line lw 2 lt 2 title "integral"\
	' | tee test.gpi | pyxplot


test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
