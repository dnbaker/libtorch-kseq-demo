

CXX?=g++

OPT+=-O3 -march=native -fopenmp
WARNING+=-Wall -Wextra -Wno-unused-parameter
LIBPATHS+= \
	-L/home-1/dbaker49@jhu.edu/miniconda3/lib/\
	-L/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch-build-ninja/lib \
	-L. -Llib64
LINKS+=-lc10 -ltorch -ltorch_cpu -lcblas -lz
LIBTORCH+=
INCLUDE+=-I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch/torch/csrc/ -Iinclude \
	-I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch/torch/csrc/api/include/ -I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch/ -I. \
	-I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch-build-ninja/aten/src/ -I/net/langmead-bigmem-ib.bluecrab.cluster/storage/dnb/code2/tch-rs/pytorch/aten/src/ -I$(LIBTORCH) \
	-Iinclude -I.
	

%: bin/%.cpp
	$(CXX) $(INCLUDE) $(WARNING) $(OPT) lib/libc10d.a $< -o $@ $(LIBPATHS) $(LINKS) -I. -Iinclude
