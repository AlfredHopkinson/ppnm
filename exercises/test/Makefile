CC = gcc #The C compiler
CFLAGS = -O -std=gnu11 # options for the C compiler
LDLIBS = -lm # libraries to link

default: out.txt    # default target is to make out.txt
	cat out.txt # show out.txt on the screen

out.txt: main            # out.txt needs hello
	./main > out.txt # run hello, output goes into out.txt

main: main.o
	cc main.o -o main
main.o: main.c
	cc -c main.c -o main.o
