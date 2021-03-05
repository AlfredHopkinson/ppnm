CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

ms.pdf: ms.tex fig-gpl.tex fig-pyxplot.pdf
	pdflatex $<
	pdflatex $<
#	bibtex $<
#	pdflatex $<

out.data.txt:main
	./$< > $@

main:main.o mygamma.o

fig-gpl.tex : fig-gpl.gpi out.data.txt
	gnuplot $<

fig-pyxplot.pdf: fig-pyxplot.gpi out.data.txt
	pyxplot $<
