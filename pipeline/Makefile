.PHONY: all help clean cleanall view report cpp clean-cpp

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

.SUFFIXES:

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

##############
## Commands ##
##############

all:
	snakemake -p -c all --rerun-incomplete --keep-going --show-failed-logs

cpp: _build_dists_v3_cpp ## Build the C++ distance builder

_build_dists_v3_cpp: _build_dists_v3.cpp ## Compile the C++ distance builder
	$(CXX) $(CXXFLAGS) $(OMP_CXXFLAGS) $(OMP_CPPFLAGS) $< -o $@ $(OMP_LDFLAGS)

help: ## Print help messages
	@echo -e "$$(grep -hE '^\S*(:.*)?##' $(MAKEFILE_LIST) \
		| sed \
			-e 's/:.*##\s*/:/' \
			-e 's/^\(.*\):\(.*\)/   \\x1b[36m\1\\x1b[m:\2/' \
			-e 's/^\([^#]\)/\1/g' \
			-e 's/: /:/g' \
			-e 's/^#\(.*\)#/\\x1b[90m\1\\x1b[m/' \
		| column -c2 -t -s : )"

report:
	snakemake --report report.html

clean: clean-cpp ## Clean
	rm -fv \
		02_order_randomization/*.txt \
		03_selection/*.txt \
		04_indices/{*.fur,*.time,*.time.tmp} \
		05_dumps/{*.fa,*.txt,*.time,*.time.tmp} \
		07_ham_dists/{*.txt,*.xz,*.time,*.time.tmp} \
		08_tsp/*.{time,time.tmp} \
		09_trees/{*.nw,*.time,*.time.tmp} \
		10_orders/{*.txt,*.time,*.time.tmp} \
		11_runs/{*.runs,*.time,*.time.tmp}
	rm -vfr 04_indices/*/ 08_tsp/*/

clean-cpp:
	rm -f _build_dists_v3_cpp

cleanall: clean ## Clean all

view: ## View output
