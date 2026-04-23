SHELL       := /bin/bash
BDM_ENV     := $(HOME)/Documents/dev/MyBDM/biodynamo/build/bin/thisbdm.sh
PYTHON      := python3
DATA_DIR    := data-export
RESULTS_DIR := results

.PHONY: all run metrics plot archive status sweep clean help

all: run metrics    ## full pipeline: simulate → metrics

run:                ## build and run simulation (~93 simulated days)
	@echo ">>> Building and running simulation..."
	source $(BDM_ENV) && biodynamo run
	@echo ">>> Done. Output: $(DATA_DIR)/populations.csv"

metrics: $(DATA_DIR)/populations.csv  ## compute fit metrics vs. scaled reference CSVs
	@echo ">>> Computing fit metrics..."
	source $(BDM_ENV) && $(PYTHON) $(DATA_DIR)/calc-error.py

plot: $(DATA_DIR)/populations.csv     ## open comparison plots (needs a display)
	source $(BDM_ENV) && $(PYTHON) $(DATA_DIR)/create-plots.py

archive:            ## save current run:  make archive NAME=baseline
	@test -n "$(NAME)" || { echo "Usage: make archive NAME=<run_name>"; exit 1; }
	mkdir -p $(RESULTS_DIR)/$(NAME)
	cp $(DATA_DIR)/populations.csv         $(RESULTS_DIR)/$(NAME)/
	cp $(DATA_DIR)/fit_metrics_summary.csv $(RESULTS_DIR)/$(NAME)/ 2>/dev/null || true
	@echo ">>> Archived to $(RESULTS_DIR)/$(NAME)/"

sweep:              ## run parameter sweep:  make sweep CONFIG=scripts/sweep_config.yaml
	@test -n "$(CONFIG)" || { echo "Usage: make sweep CONFIG=scripts/sweep_config.yaml"; exit 1; }
	source $(BDM_ENV) && $(PYTHON) scripts/run_sweep.py --config $(CONFIG)

status:             ## show current simulation state and fit metrics
	@echo ""
	@echo "=== Pancreatic Tumor ABM — Status ==="
	@echo ""
	@if [ -f $(DATA_DIR)/populations.csv ]; then \
		echo "  Last run  : $$(stat -c '%y' $(DATA_DIR)/populations.csv | cut -d. -f1)"; \
		echo "  Days sim  : $$(tail -1 $(DATA_DIR)/populations.csv | awk -F, '{printf "%s", $$2}')"; \
		echo "  Agents    : $$(tail -1 $(DATA_DIR)/populations.csv | awk -F, '{printf "%s", $$9}')"; \
	else \
		echo "  populations.csv not found — run: make run"; \
	fi
	@echo ""
	@if [ -f $(DATA_DIR)/fit_metrics_summary.csv ]; then \
		echo "  Fit metrics (vs. scaled reference):"; \
		awk -F, 'NR>1 { printf "    %-8s  R²=%6.3f  MAPE=%6.1f%%\n", $$1, $$10, $$9 }' $(DATA_DIR)/fit_metrics_summary.csv; \
	else \
		echo "  fit_metrics_summary.csv not found — run: make metrics"; \
	fi
	@echo ""
	@echo "  Archived runs:"
	@ls $(RESULTS_DIR)/ 2>/dev/null | grep -v '\.gitkeep' | sed 's/^/    /' || echo "    (none yet — use: make archive NAME=<name>)"
	@echo ""

clean:              ## remove generated outputs (keeps archived results)
	rm -f $(DATA_DIR)/populations.csv $(DATA_DIR)/fit_metrics_summary.csv

help:               ## show this help
	@echo ""
	@echo "Pancreatic Tumor ABM — available targets:"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?##' $(MAKEFILE_LIST) | \
		awk 'BEGIN{FS="##"} { split($$1,t,":"); printf "  \033[36m%-12s\033[0m%s\n", t[1], $$2}'
	@echo ""
