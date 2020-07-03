#ARCH ?= #k100#geocluster #gpupc1 #D
#MPI_ON ?= 1
#ARCH ?= Nastya

ifdef MPI_ON
GCC  ?= mpic++ 
else
ifeq ($(ARCH),Nastya)
GCC  ?= g++-7
else 
GCC  ?= g++
endif
endif

ifeq ($(ARCH),k60)
NVCC := nvcc -std=c++11 -ccbin $(GCC) -O3 -I/common/NVIDIA_CUDA-10.0_Samples/common/inc/ -I./
GENCODE_SM := -arch=sm_70
else ifeq ($(ARCH),k100)
NVCC := /common/cuda-6.5/bin/nvcc -ccbin $(GCC) -O3 
GENCODE_SM := -arch=sm_20
else ifeq ($(ARCH),D)
NVCC := /home/zakirov/cuda-7.5/bin/nvcc -ccbin $(GCC) -O3 -g 
GENCODE_SM := -arch=sm_35
else ifeq ($(ARCH),photon)
NVCC := /mnt/D/home/zakirov/cuda-7.5/bin/nvcc -ccbin $(GCC) -O3 
GENCODE_SM := -arch=sm_50
else ifeq ($(ARCH),supermic)
NVCC := /usr/local/cuda/bin/nvcc -ccbin $(GCC) -O3 -std=c++14 -g -I/usr/local/cuda/include/ -I./
GENCODE_SM := -arch=sm_61
else ifeq ($(ARCH),plasma)
NVCC := /usr/local/cuda-10.0/bin/nvcc -ccbin $(GCC) -O3 -std=c++14 --expt-relaxed-constexpr -g -I/usr/local/cuda-10.0/include/ -I./ -I./Sprout/
GENCODE_SM := -arch=sm_52
else ifeq ($(ARCH),ion)
NVCC := /mnt/D/home/zakirov/cuda-7.5/bin/nvcc -ccbin $(GCC) -O3 
GENCODE_SM := -arch=sm_50
else ifeq ($(ARCH),kiae)
NVCC := nvcc -ccbin $(GCC) -O3
GENCODE_SM := -arch=sm_37
NOG=e1
else
NVCC := nvcc -ccbin $(GCC) -O3  
GENCODE_SM := -arch=sm_61
endif 

ALL_DEFS := SET_FROM_PYTHON NX NY NZ USE_FLOAT USE_DOUBLE NOGL
CDEFS := $(foreach f, $(ALL_DEFS), $(if $($f),-D$f=$($f)))

# internal flags
NVCCFLAGS   := -Xptxas="-v"   -Xcudafe "--diag_suppress=declared_but_not_referenced" --expt-extended-lambda
CCFLAGS     := -Ofast -fopenmp -fPIC $(CDEFS) 
NVCCLDFLAGS :=
LDFLAGS     := -L./
ifeq ($(ARCH),supermic)
LDFLAGS     := -L./ -L/usr/local/cuda/lib64
endif 
ifeq ($(ARCH),kiae)
CCFLAGS     := -O3 -fopenmp -fPIC $(CDEFS) 
endif

# Extra user flags
EXTRA_NVCCFLAGS   ?=
EXTRA_NVCCLDFLAGS ?=
EXTRA_LDFLAGS     ?=
EXTRA_CCFLAGS     ?= -I./Sprout/ #-std=c++11

ifeq ($(ARCH),kiae)
INCLUDES  := -I./
LIBRARIES := -lcudart -lglut -lGL -lcufft       -lgomp -lpthread 
else
INCLUDES  := -I./ 
LIBRARIES := -lcudart -lglut -lGL -lcufft -lpng -lgomp -lpthread
endif
ifdef MPI_ON
LIBRARIES := -lmpi $(LIBRARIES)
endif

################################################################################
GENCODE_SM20  := #-gencode arch=compute_20,code=sm_21
GENCODE_SM30  := #-gencode arch=compute_30,code=sm_30
GENCODE_SM35  := #-gencode arch=compute_35,code=\"sm_35,compute_35\"
GENCODE_SM50  := #-arch=sm_20
GENCODE_FLAGS := $(GENCODE_SM50) $(GENCODE_SM35) $(GENCODE_SM30) $(GENCODE_SM20) $(GENCODE_SM) 
ALL_CCFLAGS   := --compiler-options="$(CCFLAGS) $(EXTRA_CCFLAGS)" 
ALL_LDFLAGS   := --linker-options="$(LDFLAGS) $(EXTRA_LDFLAGS)"
################################################################################
