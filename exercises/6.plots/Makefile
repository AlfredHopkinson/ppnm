CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gls-config --libs)

data.txt: main
	./$< > $@



test:
	echo $(CFLAGS)
	echo $(LDLIBS)
