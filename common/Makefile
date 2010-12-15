
OIFITS_DIR = ../lib/oifitslib
OIFITS_INC = -I$(OIFITS_DIR)
OIFITS_LIB = -L$(OIFITS_DIR)
OIFITS_OBJ = $(OIFITS_DIR)/read_fits.o $(OIFITS_DIR)/write_fits.o \
            $(OIFITS_DIR)/free_fits.o $(OIFITS_DIR)/oifile.o \
            $(OIFITS_DIR)/oifilter.o $(OIFITS_DIR)/oicheck.o \
            $(OIFITS_DIR)/oimerge.o

# Create a list of all of the .cpp files:
CPP_SOURCES = $(wildcard *.cpp)
CPP_OBJECTS = $(patsubst %.cpp, %.o, $(CPP_SOURCES))

# Create a list of all of the .c files:
CC_SOURCES = $(wildcard *.c)
CC_OBJECTS = $(patsubst %.c, %.o, $(CC_SOURCES))

GTK_INC = `pkg-config --cflags gtk+-2.0`
GTK_LIB = `pkg-config --libs gtk+-2.0`

INCLUDES = $(OIFITS_INC) $(GTK_INC)
LIBS = $(OIFITS_LIB) $(GTK_LIB)
LDFLAGS = -loifits 

# Compilers:
CPP = g++
CPPFLAGS = -ggdb -g -Wall --pedantic -Wextra $(INCLUDES)

CC = gcc
CCFLAGS = -g -Wall


############# Make Targets below #############

default: cc_common cpp_common libOISim.a

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $<
	
%.o : %.c
	$(CC) $(CCFLAGS) -c $<
	
cc_common: $(CC_OBJECTS)

cpp_common: $(CPP_OBJECTS)
	
libOISim.a:	cc_common cpp_common
	ar -ruc $@ $(OIFITS_OBJ) $(CPP_OBJECTS) $(CC_OBJECTS)

	
# CLEANER
.PHONY: clean
clean:
	rm -f *.o
	rm -f *.a
