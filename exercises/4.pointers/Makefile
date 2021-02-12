CC =gcc
CFLAGS = -std=gnu11  -Wall
LDFLAGS = -I.
LDLIBS = -lm

out.txt: main
	./$< > $@  #$< is the first variable so main target is another varibale which is out.txt
	.vector >> $@

clean:
	$(RM) main out.txt
