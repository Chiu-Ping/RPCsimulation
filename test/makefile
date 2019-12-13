OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR) -I/usr/include/c++/5

# Debug flags
#CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm
#LDFLAGS += -g

rpc1: rpc1.c 
	$(CXX) $(CFLAGS) rpc1.c
	$(CXX) -o rpc1 rpc1.o $(LDFLAGS)
	rm rpc1.o
	
pi: pi.c 
	$(CXX) $(CFLAGS) pi.c
	$(CXX) -o pi pi.o $(LDFLAGS)
	rm pi.o

gasfile: gasfile.C 
	$(CXX) $(CFLAGS) gasfile.C
	$(CXX) -o gasfile gasfile.o $(LDFLAGS)
	rm gasfile.o