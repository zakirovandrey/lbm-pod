include header.mk

export ALL_DEFS

# Target rules
all: build 

build: _lbm.so _cuda_pyhelper.so swig

swig: _cuda_pyhelper.so cuda_pyhelper.i cuda_pyhelper.h
cuda_pyhelper_wrap.cu: cuda_pyhelper.i cuda_pyhelper.h 
	swig -python -globals G -c++ -I/usr/local/cuda/include/ -o cuda_pyhelper_wrap.cu $<
cuda_pyhelper_wrap.o: cuda_pyhelper_wrap.cu cuda_pyhelper.h
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -dc $< -I/usr/include/python/ -I/usr/include/python2.7/ -I/usr/include/python2.6/ -o $@ \
          -Xcudafe "--diag_suppress=set_but_not_used" -w
_cuda_pyhelper.so: cuda_pyhelper_wrap.o
	$(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) $(LDFLAGS) -L./ $+ -o $@ $(LIBRARIES) -shared

im3D.o: im3D.cu im3D.hpp cuda_math.h fpal.h im2D.h Arr3Dpars.hpp ArrLst.hpp 
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<
err.o: err.cu im3D.hpp Arr3Dpars.hpp err.h
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<

DATA_FILE = data.cu 
RUN_FILE = run.cu 
UPDATE_FILE = update.cu
DROP_FILE = drop.cu
INIT_FILE = init.cu
DRAW_FILE = draw_im3D.cu

ker-run.dep      ker-run.o  : $(DATA_FILE)
data.dep         data.o     : $(RUN_FILE)
step-update.dep  update.o   : $(UPDATE_FILE)
drop.dep         drop.o     : $(DROP_FILE)
init.dep         init.o     : $(INIT_FILE)
draw_im3D.dep    draw_im3D.o: $(DRAW_FILE)

ker-run.dep data.dep step-update.dep drop.dep init.dep draw.dep:
	@echo "\033[0;32m   ===== Building dependencies file for $(current_dir)/$^  \033[0m"
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ >$@ ;

ifneq ($(filter clean,$(MAKECMDGOALS)),clean)
include ker-run.dep data.dep step-update.dep drop.dep init.dep draw_im3D.dep
endif

ker-run.o data.o update.o drop.o init.o draw_im3D.o:
	@echo "\033[0;32m  |: Compiling $(current_dir)/$<  \033[0m"
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<

common.o: ker-run.o data.o update.o drop.o init.o draw_im3D.o
	ar rvs $@ $+

lbm_wrap.cxx: lbm.i phys.h
	swig -python -globals G -c++ -o lbm_wrap.cxx $<
lbm_wrap.o: lbm_wrap.cxx phys.h
	$(GCC) $(INCLUDES) $(CCFLAGS) -c $< -fPIC -I/usr/include/python/ -I/usr/include/python2.7/ -I/usr/include/python2.6/ -DIGNORE_CUDA_EXTENSIONS=1 -o $@

all_lbm.so: common.o im3D.o err.o
	@echo "\033[0;32m  |::::  Linking the target $@ \033[0m"
	$(EXEC) $(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) $(LDFLAGS) -o $@ $(filter %.o,$+) $(LIBRARIES) --shared
	@echo "  ---  Linking finished --- "

_lbm.so: all_lbm.so lbm_wrap.o
	@echo "\033[0;32m  |::::  Linking the target $@ \033[0m"
	$(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) $(LDFLAGS) -L./ $< lbm_wrap.o -o $@ $(LIBRARIES) -shared
	@echo "  ---  Linking finished --- "

#lbm: ker_fastlab.o fill_init.o update.o moments.o drop.o mesh_reform.o draw_im3D.o ../core/core_fastlab.o
#	$(EXEC) $(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) $(LDFLAGS) -o $@ $+ $(LIBRARIES)

run: lbm
	$(EXEC) ./lbm

clean:
	$(EXEC) rm -f lbm_wrap.o all_lbm.so _lbm.so cuda_pyhelper_wrap.o _cuda_pyhelper.so ../cuda_pyhelper_wrap.cu \
	*.dep \
	*.o lbm lbm.py* ./*_wrap.cxx
clobber: clean

.PHONY: clean all

#Xcudafe --diag_suppress warnings list
#
