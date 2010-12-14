# MAKEFILE FOR THE MROI VISIBILITY SIMULATOR

##### May need to edit this section #####

OIFITS_DIR = ./lib/oifitslib
COMMON_DIR = ./common
VSIM_DIR = ./vis_sim

.PHONY: all
all: oifitslib common vsim

# We need to link with the OIFITS library, make it.
.PHONY: oifitslib
oifitslib:
	$(MAKE) liboifits.a -C $(OIFITS_DIR)

.PHONY: common
common:
	$(MAKE) default -C $(COMMON_DIR)
	
.PHONY: vsim
vsim:
	$(MAKE) vsim -C $(VSIM_DIR)

# CLEANER
.PHONY: clean
clean:
	rm -f *.o
	rm -f ../bin/*
	$(MAKE) clean -C $(COMMON_DIR)
	$(MAKE) clean -C $(OIFITS_DIR)
