CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt 

main: main.o ann.o mini.o

out.txt:main
	./$< > out.txt


.PHONEY:clean
clean:	
	$(RM) main.o main out.txt num.txt mini.o mini ann.o ann

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
