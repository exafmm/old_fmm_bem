.SUFFIXES: .cu .o

#NVCC = nvcc --ptxas-options=-v -use_fast_math -I$(SDK_INSTALL_PATH)/common/inc
NVCC = nvcc --ptxas-options=-v --host-compilation 'C++' -Xcompiler "-march=core2" -Xcompiler "-m64" -Xcompiler "-O3" -Xcompiler "-ffast-math" -Xcompiler "-fopenmp" --machine 64 -DUNIX -O3 -I. -I$(CUDA_INSTALL_PATH)/include -I$(SDK_INSTALL_PATH)/common/inc

.cu.o:
	$(NVCC) -c $<

libgputree.a: gputree.o
	ar ruv $@ $?
	ranlib $@
	$(RM) *.o
libgpufmm.a: gpufmm.o
	ar ruv $@ $?
	ranlib $@
	$(RM) *.o
clean:
	$(RM) libgputree.a libgpufmm.a
tree:
	make libgputree.a
fmm:
	make libgpufmm.a
all:
	make clean
	make libgputree.a
	make libgpufmm.a
