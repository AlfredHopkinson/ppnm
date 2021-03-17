CC = gcc
LDLIBS = -lm
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
 

cubic.out.txt: cubic
		./$< > $@
