.SUFFIXES: .cxx .o

KERNEL = ../kernel
FMM = ../fmm
SORT = ../sort
MPI = ../mpi
# Tree structure common to Serial & Parallel
OBJ += $(FMM)/fmm.o $(FMM)/alloc.o $(FMM)/boxdatai.o $(FMM)/boxdataj.o\
	$(FMM)/boxm2mi.o $(FMM)/boxm2mj.o\
	$(FMM)/boxn.o $(FMM)/boxn1.o $(FMM)/boxc.o\
	$(FMM)/ijbox.o $(FMM)/memory.o $(SORT)/sortvar.o $(SORT)/csort.o
# Tree structure only Serial
OBJ_SERIAL += $(FMM)/nlevel.o $(FMM)/boxallocate.o $(FMM)/boxpart.o $(FMM)/ifbox.o\
	$(SORT)/sorti.o $(SORT)/sortj.o
# Tree structure only Parallel
OBJ_PARALLEL = $(MPI)/mpisendrecv.o $(MPI)/mpialltoallv.o $(MPI)/mpishift.o\
	$(MPI)/mpinlevel.o $(MPI)/mpiboxallocate.o\
	$(MPI)/mpisendp2p.o $(MPI)/mpisendm2l.o $(MPI)/mpiboxpart.o\
	$(MPI)/mpiijcnt.o $(MPI)/mpiifbox.o $(MPI)/mpirange.o\
	$(SORT)/mpisorti.o $(SORT)/mpisortj.o\
	$(SORT)/mpisortvar.o

# Kernels common to CPU & GPU
OBJ += $(KERNEL)/precalc.o $(KERNEL)/spharot.o $(KERNEL)/multipole.o\
	$(KERNEL)/multipoled.o $(KERNEL)/factorial.o
# Kernels only CPU
OBJ_CPU += $(KERNEL)/G_p2p.o $(KERNEL)/Gni_p2p.o $(KERNEL)/Gnj_p2p.o\
	$(KERNEL)/G_p2m.o $(KERNEL)/Gn_p2m.o\
	$(KERNEL)/m2m.o $(KERNEL)/m2l.o $(KERNEL)/l2l.o\
	$(KERNEL)/G_l2p.o $(KERNEL)/Gn_l2p.o\
	$(KERNEL)/G_dir.o $(KERNEL)/Gni_dir.o $(KERNEL)/Gnj_dir.o
# Kernels only GPU
OBJ_GPU += $(KERNEL)/G_p2p_gpu.o $(KERNEL)/Gni_p2p_gpu.o $(KERNEL)/Gnj_p2p_gpu.o\
	$(KERNEL)/G_p2m_gpu.o $(KERNEL)/Gn_p2m_gpu.o\
	$(KERNEL)/m2m_gpu.o $(KERNEL)/m2l_gpu.o $(KERNEL)/l2l_gpu.o\
	$(KERNEL)/G_l2p_gpu.o $(KERNEL)/Gn_l2p_gpu.o\
	$(KERNEL)/G_dir_gpu.o $(KERNEL)/Gni_dir_gpu.o $(KERNEL)/Gnj_dir_gpu.o
# BEM routines
OBJ += geometry.o bem3d.o
# GPU libraries
GPU_LIB = ../libgpu/libgpufmm.a -L$(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart

# Make commands
cpu: main.o $(OBJ) $(OBJ_SERIAL) $(OBJ_CPU)
	$(CLINKER) $? $(PETSC_SNES_LIB)
	@$(RM) *.o
gpu: main.o $(OBJ) $(OBJ_SERIAL) $(OBJ_GPU)
	$(CLINKER) $? $(PETSC_SNES_LIB) $(GPU_LIB)
	@$(RM) *.o
cpumpi: mpimain.o $(OBJ) $(OBJ_PARALLEL) $(OBJ_CPU)
	$(CLINKER) $? $(PETSC_SNES_LIB)
	@$(RM) *.o
gpumpi: mpimain.o $(OBJ) $(OBJ_PARALLEL) $(OBJ_GPU)
	$(CLINKER) $? $(PETSC_SNES_LIB) $(GPU_LIB)
	@$(RM) *.o

cleanall:
	@$(RM) $(DEF) $(PAR) *.o *.out *.dat $(KERNEL)/*.o $(FMM)/*.o $(SORT)/*.o $(MPI)/*.o

save:
	tar zcvf ../../bem.tgz ../../bem
ifeq (0,0)
# Use this one
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
else
# Rio sometimes likes to use intel compilers for debugging
CLINKER = mpicxx -check -traceback -g -I/nfs/home/yokota/petsc-dev/linux-gnu-cxx-opt/include -I/nfs/home/yokota/petsc-dev/include -I/nfs/home/yokota/petsc-dev/linux-gnu-cxx-opt/include -I/nfs/home/yokota/petsc-dev/include/sieve -I/usr/include
CLINKER = mpicxx -O3 -I$(PETSC_DIR)/linux-gnu-cxx-opt/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/linux-gnu-cxx-opt/include -I$(PETSC_DIR)/include/sieve -I/usr/include
PETSC_SNES_LIB = -Wl,-rpath,$(PETSC_DIR)/linux-gnu-cxx-opt/lib -L$(PETSC_DIR)/linux-gnu-cxx-opt/lib -lpetsc  -lX11 -Wl,-rpath,$(PETSC_DIR)/linux-gnu-cxx-opt/lib -L$(PETSC_DIR)/linux-gnu-cxx-opt/lib -lparmetis -lmetis -lflapack -lfblas -lnsl -lrt -Wl,-rpath,/usr/local/cuda/lib -L/usr/local/cuda/lib -Wl,-rpath,/nfs/home/yokota/fftw/lib -L/nfs/home/yokota/fftw/lib -Wl,-rpath,$(PETSC_DIR) -L$(PETSC_DIR) -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -ldl -lmpich -lpthread -lrt -lgcc_s -lmpichf90 -lgfortran -lm -lmpichcxx -lstdc++
.cxx.o:
	$(CLINKER) -c $< -o $@
endif
