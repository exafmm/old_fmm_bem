.SUFFIXES: .cxx

CC = gcc -O3 -Wall -ffast-math -funroll-loops -march=nocona

.c.o:
	$(CC) -c $<

libcpufmm.a: cpufmm.o
	ar ruv $@ $?
	ranlib $@
	$(RM) *.o
clean:
	$(RM) libcpufmm.a
all:
	make clean
	make libcpufmm.a
run:
	make all
	cd ../test;make fcpu;./a.out
