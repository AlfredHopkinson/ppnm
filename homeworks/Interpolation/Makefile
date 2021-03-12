CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


main: main.c lin.c bin_search.c

lin.out.txt:lin
	./$< > lin.out.txt

.PHONEY:clean
clean:
	$(RM) main.o main out.txt

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
