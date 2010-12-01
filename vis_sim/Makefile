# MAKEFILE FOR THE MROI VISIBILITY SIMULATOR

##### May need to edit this section #####

# Library Paths:
CFITSIO_INC = -I/usr/include -I/usr/local/include
CFITSIO_LIB = -L/usr/lib -L/usr/local/lib -lcfitsio

OIFITS_DIR = ../lib/oifitslib
OIFITS_INC = -I$(OIFITS_DIR)
OIFITS_LIB = -L$(OIFITS_DIR)

COMMON_DIR = ../common
COMMON_INC = -I$(COMMON_DIR)
COMMON_LIB = -L$(COMMON_DIR)

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))

INCLUDES = $(COMMON_INC) $(CFITSIO_INC) $(OIFITS_INC) \
	`pkg-config --cflags fftw3` `pkg-config --cflags gtk+-2.0`
LDFLAGS = `pkg-config --libs fftw3` `pkg-config --libs gtk+-2.0` -loifits -lOISim
LIBS = $(CFITSIO_LIB) $(OIFITS_LIB) $(COMMON_LIB)

# Compilers/utilities to use
CPP = g++ 

# Compiler flags
CPPFLAGS = -ggdb -g -Wall --pedantic -Wextra $(INCLUDES)

EXECUTABLE = oifits-sim

##### end of user-editable section #####

# Rule to make .o file from .cpp file
%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $<

.PHONY: all
all: oifitslib common vsim

# We need to link with the OIFITS library, make it.
.PHONY: oifitslib
oifitslib:
	$(MAKE) liboifits.a -C $(OIFITS_DIR)

.PHONY: common
common:
	$(MAKE) default -C $(COMMON_DIR)
	
%.o : %.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@
	
vsim: $(OBJECTS)
	$(CPP) -o $(EXECUTABLE) $(OBJECTS) $(CPPFLAGS) $(LIBS) $(LDFLAGS)

# CLEANER
.PHONY: clean
clean:
	rm -f *.o
	rm -f ../bin/*
	$(MAKE) clean -C $(COMMON_DIR)
	$(MAKE) clean -C $(OIFITS_DIR)
