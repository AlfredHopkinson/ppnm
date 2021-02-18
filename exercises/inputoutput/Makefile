CC = gcc
CFLAGS =
LDLIBS = -lm

out.txt: inoutputA inoutputB inoutputC
		
inoutputA.o: inoutputA.c
		$(CC)$(CFLAGS) -c inoutputA.c -o inoutputA.o
inoutputB.o: inoutputB.c
		$(CC)$(CFLAGS) -c inoutputB.c -o inoutputB.o
inoutputC.o: inoutputC.c
	                $(CC)$(CFLAGS) -c inoutputC.c -o inoutputC.o

inoutputA: inoutputA.o
		$(CC) $(LDFLAGS) inoutputA.o -o inoutputA $(LDLIBS)
inoutputB: inoutputB.o
	                $(CC) $(LDFLAGS) inoutputB.o -o inoutputB $(LDLIBS)
inoutputC: inoutputC.o
	                $(CC) $(LDFLAGS) inoutputC.o -o inoutputC $(LDLIBS)


.PHONEY:clean
clean:
