# ================================================
# Use this makefile to have frequent used command
# and use them through "make 'command'"
# Author: Tony Liang
# ================================================

# Help message variable
define HELP_MESSAGE
This is the default help message to display.

Usage: make [target]

Target:
  - help          Display this help message (default prints it if not target provided)
  - clean         Delete a lots of stuff like log files and work/ output

endef
export HELP_MESSAGE

# By default should run this help target
help:
	@echo "$${HELP_MESSAGE}"

# Use this target to clean intermediate files
.PHONY: clean
clean:
	@echo "Removing intermediate files like logs and work/ ..."
	@rm -rf work/
	@rm -f .nextflow.log*
	@rm -rf .nextflow/
	@rm -rf results/