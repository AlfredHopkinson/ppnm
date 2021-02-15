CC =gcc
CFLAGS = 
LDLIBS = -lm

out.txt: epsilon
	./$< > $@  #$< is the first variable so main target is another varibale which is out.txt

epsilon.o: epsilon.c
	$(CC)$(CFLAGS) -c epsilon.c -o epsilon.o
equal.o: equal.c
	$(CC)$(CFLAGS) -c equal.c -o equal.o

epsilon: epsilon.o equal.o
	$(CC) $(LDFLAGS) epsilon.o equal.o -o epsilon $(LDLIBS)

clean:
	$(RM) epsilon out.txt epsilon.o equal.o
