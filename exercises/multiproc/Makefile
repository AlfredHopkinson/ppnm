CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp


default: out.main.txt out.omp.txt

out.main.txt:main
		./$< > $@
out.omp.txt:openmp
		./$< > $@
