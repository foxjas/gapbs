# See LICENSE.txt for license details.

# CXX_FLAGS += -g -std=c++11 -O3 -Wall -march=native -mavx512f



CXX_FLAGS += -g -std=c++11 -O3 -Wall -march=native -xMIC-AVX512 -qopt-report=5
CXX = icpc

PAR_FLAG = -fopenmp

ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG = -qopenmp
endif

ifneq (,$(findstring sunCC,$(CXX)))
	CXX_FLAGS = -std=c++11 -xO3 -m64 -xtarget=native
	PAR_FLAG = -xopenmp
endif

ifneq ($(SERIAL), 1)
	CXX_FLAGS += $(PAR_FLAG)
endif

#KERNELS = bc bfs cc pr sssp tc bc_lrb pr_lrb
KERNELS = pr pr_lrb 

SUITE = $(KERNELS) 
#SUITE = $(KERNELS) converter

.PHONY: all
all: $(SUITE)

% : src/%.cc src/*.h
	$(CXX) $(CXX_FLAGS) $< -o $@

# Testing
include test/test.mk

# Benchmark Automation
include benchmark/bench.mk


.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*
