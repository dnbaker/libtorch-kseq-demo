

CXX?=g++
NVCC?=nvcc

OPT+=-O3 #-fopenmp
WARNING+= #-Wall -Wextra -Wno-unused-parameter
LIBPATHS+= \
	-L~/miniconda3/envs/cuda/lib -L/home-1/dbaker49@jhu.edu/miniconda3/envs/cuda/lib/python3.9/site-packages/torch/lib/ #-L. -L64
LINKS+=-lc10 -ltorch -ltorch_cpu -lcblas -lz
LIBTORCH+=
INCLUDE+=-I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch/torch/csrc/ -Iinclude \
	-I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch/torch/csrc/api/include/ -I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch/ -I. \
	-I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch-build-ninja/aten/src/ -I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch/aten/src/ -I$(LIBTORCH) \
	-Iinclude -I.
	

%-gpu: bin/%.cpp
	$(NVCC) $(INCLUDE) $(WARNING) $(OPT) lib/libc10d.a $< -o $@ $(LIBPATHS) $(LINKS) -I. -Iinclude -ltorch_cuda -lc10_cuda
