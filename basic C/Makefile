C = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm # Math Libary


out.txt: hello math
	./hello > out.txt
	./math >> out.txt # Double >> means append

hello.o: hello.c hello.h
	$(CC) $(CFLAGS) -c hello.c -o hello.o

hello: hello.o world.o 
	$(CC) $(LDFLAGS) hello.o world.o -o hello $(LDLIBS)

world.o: world.c hello.h
	$(CC) $(CLFAGS) -c $< -o $@

math: math.o

clean:
	$(RM) *.o out* hello

