NVCC = /opt/cuda/bin/nvcc   
LIBS = -lcufft -lcublas -lhdf5 -lhdf5_hl -lcurand -lcusparse
DEBUG = -g
GPU_SOURCES = $(wildcard src/*.cu)
GPU_OBJECTS = $(GPU_SOURCES:.cu=.o)



all: $(GPU_OBJECTS)
	$(NVCC) -o taylorC $(GPU_OBJECTS) $(LIBS)

$(GPU_OBJECTS): src/%.o: src/%.cu
	$(NVCC) -c   $< -o $@

clean:
	rm src/*.o taylorC
