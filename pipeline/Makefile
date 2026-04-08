.PHONY: all help clean cleanall view report

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

.SUFFIXES:

##############
## Commands ##
##############

all:
	snakemake -p -c all --rerun-incomplete --keep-going --show-failed-logs

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

clean: ## Clean
	rm -fv \
		02_order_randomization/*.txt \
		03_selection/*.txt \
		04_indices/{*.fur,*.time,*.time.tmp} \
		05_dumps/{*.fa,*.txt,*.time,*.time.tmp} \
		07_ham_dists/{*.txt,*.xz,*.time,*.time.tmp} \
		08_tsp/*.{time,time.tmp} \
		08_tsp/*/*.{time,time.tmp} \
		09_trees/{*.nw,*.time,*.time.tmp} \
		10_orders/{*.txt,*.time,*.time.tmp} \
		11_runs/{*.runs,*.time,*.time.tmp}
	rm -vfr 04_indices/*/ 08_tsp/*/

cleanall: clean ## Clean all

view: ## View output
