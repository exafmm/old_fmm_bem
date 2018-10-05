.SUFFIXES: .cxx .o

FMM = ../fmm
KERNEL = ../kernel
MPI = ../mpi
SORT = ../sort
TREE = ../tree
# FMM tree structure common to Serial & Parallel
OBJ_FMM += $(FMM)/fmm.o $(FMM)/alloc.o $(FMM)/boxdatai.o $(FMM)/boxdataj.o\
	 $(FMM)/boxm2mi.o $(FMM)/boxm2mj.o\
	$(FMM)/boxn.o $(FMM)/boxn1.o $(FMM)/boxc.o\
	$(FMM)/ijbox.o $(FMM)/memory.o $(SORT)/sortvar.o $(SORT)/csort.o
# FMM tree structure only Serial
OBJ_FMM_SERIAL += $(FMM)/nlevel.o $(FMM)/boxallocate.o $(FMM)/boxpart.o $(FMM)/ifbox.o\
	$(SORT)/sorti.o $(SORT)/sortj.o
# FMM tree structure only Parallel
OBJ_FMM_PARALLEL = $(MPI)/mpisendrecv.o $(MPI)/mpialltoallv.o $(MPI)/mpishift.o\
	$(MPI)/mpinlevel.o $(MPI)/mpiboxallocate.o\
	$(MPI)/mpisendp2p.o $(MPI)/mpisendm2l.o $(MPI)/mpiboxpart.o\
	$(MPI)/mpiijcnt.o $(MPI)/mpiifbox.o $(MPI)/mpirange.o\
	$(SORT)/mpisorti.o $(SORT)/mpisortj.o\
	$(SORT)/mpisortvar.o

# FMM kernels common to CPU & GPU
OBJ_FMM += $(KERNEL)/precalc.o $(KERNEL)/precalch.o $(KERNEL)/spharot.o $(KERNEL)/multipole.o\
	$(KERNEL)/multipoled.o $(KERNEL)/multipoleh.o $(KERNEL)/factorial.o $(KERNEL)/in.o
# FMM kernels only CPU
OBJ_FMM_CPU += $(KERNEL)/G_p2p_cpu.o $(KERNEL)/Gni_p2p.o $(KERNEL)/Gnj_p2p.o\
	$(KERNEL)/G_p2m.o $(KERNEL)/Gn_p2m.o\
	$(KERNEL)/m2m.o $(KERNEL)/m2l.o $(KERNEL)/l2l.o\
	$(KERNEL)/G_l2p.o $(KERNEL)/Gn_l2p.o\
	$(KERNEL)/G_dir_cpu.o $(KERNEL)/Gni_dir.o $(KERNEL)/Gnj_dir.o\
	$(KERNEL)/Helmholtz_G_p2p.o $(KERNEL)/Helmholtz_G_p2m.o $(KERNEL)/Helmholtz_m2m.o\
	$(KERNEL)/Helmholtz_m2l.o $(KERNEL)/Helmholtz_G_l2p.o\
	$(KERNEL)/Helmholtz_G_dir.o
# FMM kernels only GPU
OBJ_FMM_GPU += $(KERNEL)/G_p2p_gpu.o $(KERNEL)/Gni_p2p_gpu.o $(KERNEL)/Gnj_p2p_gpu.o\
	$(KERNEL)/G_p2m_gpu.o $(KERNEL)/Gn_p2m_gpu.o\
	$(KERNEL)/m2m_gpu.o $(KERNEL)/m2l_gpu.o $(KERNEL)/l2l_gpu.o\
	$(KERNEL)/G_l2p_gpu.o $(KERNEL)/Gn_l2p_gpu.o\
	$(KERNEL)/G_dir_gpu.o $(KERNEL)/Gni_dir_gpu.o $(KERNEL)/Gnj_dir_gpu.o

# Treecode tree structure common to Serial & Parallel
OBJ_TREE = $(TREE)/tree.o $(TREE)/distribute.o $(TREE)/time.o $(FMM)/memory.o $(FMM)/alloc.o\
	$(KERNEL)/G_dir_gpu.o $(KERNEL)/Gni_dir.o $(KERNEL)/Gnj_dir.o
# Treecode tree structure only Serial
OBJ_TREE_SERIAL = $(TREE)/nbody_MP.o
# Treecode tree structure only Parallel
OBJ_TREE_PARALLEL = $(TREE)/nbody_MP_MPI.o $(MPI)/mpirange.o $(MPI)/mpishift.o

# Treecode kernels only CPU
OBJ_TREE_CPU = $(TREE)/gravity.o $(TREE)/BHtree.o
# Treecode kernels only GPU
OBJ_TREE_GPU = $(TREE)/gravity_g80.o $(TREE)/BHtree_g80.o

# GPU libraries
GPU_LIB = -L$(CUDA_INSTALL_PATH)/lib64 -L$(SDK_INSTALL_PATH)/lib -lcudart -lGL -lGLU -lcutil
TREE_CPU_LIB = -L../libcpu -lcputree
TREE_GPU_LIB = -L../libgpu -lgputree
FMM_CPU_LIB = -L../libcpu -lcpufmm
FMM_GPU_LIB = -L../libgpu -lgpufmm

# Make commands
tcpu: test.o $(OBJ_TREE) $(OBJ_TREE_SERIAL) $(OBJ_TREE_CPU)
	$(CLINKER) -DPHANTOM $? $(TREE_CPU_LIB) $(FMM_GPU_LIB) $(GPU_LIB)
tgpu: test.o $(OBJ_TREE) $(OBJ_TREE_SERIAL) $(OBJ_TREE_GPU)
	$(CLINKER) $? $(TREE_GPU_LIB) $(FMM_GPU_LIB) $(GPU_LIB)
tcpumpi: mpitest.o $(OBJ_TREE) $(OBJ_TREE_PARALLEL) $(OBJ_TREE_CPU)
	$(CLINKER) -DPHANTOM $? $(TREE_CPU_LIB) $(FMM_GPU_LIB) $(GPU_LIB)
tgpumpi: mpitest.o $(OBJ_TREE) $(OBJ_TREE_PARALLEL) $(OBJ_TREE_GPU)
	$(CLINKER) $? $(TREE_GPU_LIB) $(FMM_GPU_LIB) $(GPU_LIB)

fcpu: test_fmm.o $(OBJ_FMM) $(OBJ_FMM_SERIAL) $(OBJ_FMM_CPU)
	$(CLINKER) -DFMM $? $(PETSC_SNES_LIB) $(FMM_CPU_LIB)
fgpu: test_fmm.o $(OBJ_FMM) $(OBJ_FMM_SERIAL) $(OBJ_FMM_GPU)
	$(CLINKER) -DFMM $? $(PETSC_SNES_LIB) $(FMM_GPU_LIB) $(GPU_LIB)
fcpumpi: mpitest_fmm.o $(OBJ_FMM) $(OBJ_FMM_PARALLEL) $(OBJ_FMM_CPU)
	$(CLINKER) -DFMM $? $(PETSC_SNES_LIB) $(FMM_CPU_LIB)
fgpumpi: mpitest_fmm.o $(OBJ_FMM) $(OBJ_FMM_PARALLEL) $(OBJ_FMM_GPU)
	$(CLINKER) -DFMM $? $(PETSC_SNES_LIB) $(FMM_GPU_LIB) $(GPU_LIB)

test_fmm.o: test.cxx
	$(CLINKER) -DFMM -c -o test_fmm.o test.cxx
mpitest_fmm.o: mpitest.cxx
	$(CLINKER) -DFMM -c -o mpitest_fmm.o mpitest.cxx
$(TREE)/gravity.o: $(TREE)/gravity.cxx
	$(CLINKER) -DPHANTOM -c -o $(TREE)/gravity.o $(TREE)/gravity.cxx
$(TREE)/gravity_g80.o: $(TREE)/gravity.cxx
	$(CLINKER) -c -o $(TREE)/gravity_g80.o $(TREE)/gravity.cxx
$(TREE)/BHtree.o: $(TREE)/BHtree.cxx
	$(CLINKER) -DPHANTOM -c -o $(TREE)/BHtree.o $(TREE)/BHtree.cxx 
$(TREE)/BHtree_g80.o: $(TREE)/BHtree.cxx
	$(CLINKER) -c -o $(TREE)/BHtree_g80.o $(TREE)/BHtree.cxx

cleanall:
	@$(RM) $(DEF) $(PAR) *.o *.out *.dat $(KERNEL)/*.o $(TREE)/*.o $(FMM)/*.o $(SORT)/*.o $(MPI)/*.o

ifeq (0,1)
# Use this one
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
else
# Rio sometimes likes to use intel compilers for debugging
# CLINKER = mpicxx -check -traceback -g
CLINKER = mpicxx -I$(CUDA_INSTALL_PATH)/include -I$(SDK_INSTALL_PATH)/common/inc $(PGINCLUDE) -O3 -Wall -ffast-math -funroll-loops -fopenmp
PETSC_SNES_LIB = -lstdc++
.cxx.o:
	$(CLINKER) -c $< -o $@
endif
