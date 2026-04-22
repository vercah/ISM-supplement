.PHONY: all clean

ifeq ($(origin CXX), default)
CXX = g++
endif

ifeq ($(origin CXXFLAGS), undefined)
CXXFLAGS = -O3 -DNDEBUG -std=c++17 -march=native
endif

OMP_CXXFLAGS ?= -fopenmp
OMP_CPPFLAGS ?=
OMP_LDFLAGS ?=

ifneq (,$(findstring clang,$(shell $(CXX) --version 2>/dev/null | head -n 1)))
OMP_CXXFLAGS = -Xpreprocessor -fopenmp
ifneq ($(wildcard /opt/homebrew/opt/libomp/include/omp.h),)
OMP_CPPFLAGS = -I/opt/homebrew/opt/libomp/include
OMP_LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp
else ifneq ($(wildcard /usr/local/opt/libomp/include/omp.h),)
OMP_CPPFLAGS = -I/usr/local/opt/libomp/include
OMP_LDFLAGS = -L/usr/local/opt/libomp/lib -lomp
endif
endif

all: _build_dists_v5

_build_dists_v5: main.cpp
	$(CXX) $(CXXFLAGS) $(OMP_CXXFLAGS) $(OMP_CPPFLAGS) $< -o $@ $(OMP_LDFLAGS)

clean:
	rm -f _build_dists_v5
