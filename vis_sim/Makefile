# MAKEFILE FOR THE MROI VISIBILITY SIMULATOR

##### May need to edit this section #####

# Library Paths:
CFITSIO_INC = -I/usr/include -I/usr/local/include
CFITSIO_LIB = -L/usr/lib -L/usr/local/lib -lcfitsio

# The OIFITS library is merged with OISim so we only need the includes
OIFITS_DIR = ../lib/oifitslib
OIFITS_INC = -I$(OIFITS_DIR)

COMMON_DIR = ../common
COMMON_INC = -I$(COMMON_DIR)
COMMON_LIB = -L$(COMMON_DIR)

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))

PKGCONFIG = pkg-config
GTK_INC = `$(PKGCONFIG) --cflags  gtk+-2.0`
GTK_LIB = `$(PKGCONFIG) --libs  gtk+-2.0`

FFT_INC = `$(PKGCONFIG) --cflags fftw3` 
FFT_LIB = `$(PKGCONFIG) --libs fftw3`

INCLUDES = $(COMMON_INC) $(CFITSIO_INC) $(OIFITS_INC) $(FFT_INC) $(GTK_INC)
LDFLAGS = -lOISim
LIBS = $(CFITSIO_LIB) $(COMMON_LIB) $(FFT_LIB) $(GTK_LIB)

# Compilers/utilities to use
CPP = g++ 

# Compiler flags
CPPFLAGS = -ggdb -g -Wall --pedantic -Wextra $(INCLUDES)

EXECUTABLE = ../oifits-sim

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
