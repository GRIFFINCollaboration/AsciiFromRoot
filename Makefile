.EXPORT_ALL_VARIABLES:

.SUFFIXES:

.PHONY: clean all

# := is only evaluated once
SHELL 		= /bin/sh

NAME		= AsciiFromRoot

LIB_DIR 	= $(HOME)/lib
BIN_DIR		= $(HOME)/bin

ROOTLIBS     	:= $(shell root-config --libs)
ROOTINC      	:= -I$(shell root-config --incdir)

COMMON_DIR 	= ../CommandLineInterface

INCLUDES        = -I$(COMMON_DIR) -I.

LIBRARIES	= CommandLineInterface 

CC		= gcc
CXX             = g++
CPPFLAGS 	= $(ROOTINC) $(INCLUDES) -fPIC
CXXFLAGS	= -std=gnu++0x -pedantic -Wall -Wno-long-long -g -O3

LDFLAGS		= -g -fPIC

LDLIBS 		= -L$(LIB_DIR) -Wl,-rpath,/opt/gcc/lib64 $(ROOTLIBS) $(addprefix -l,$(LIBRARIES))

LOADLIBES = \
	RootUtilities.o 

# -------------------- implicit rules --------------------
# n.o made from n.c by 		$(CC) -c $(CPPFLAGS) $(CFLAGS)
# n.o made from n.cc by 	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS)
# n made from n.o by 		$(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)

# -------------------- rules --------------------

all:  AsciiFromRoot
	@echo Done

# -------------------- libraries --------------------

$(LIB_DIR)/libCommandLineInterface.so:
	@cd $(COMMON_DIR); make $@

# -------------------- pattern rules --------------------
# this rule sets the name of the .cc file at the beginning of the line (easier to find)

%.o: %.cc %.hh
	$(CXX) $< -c $(CPPFLAGS) $(CXXFLAGS) -o $@

# -------------------- default rule for executables --------------------

%: %.cc $(LOADLIBES)
	$(CXX) $< $(CXXFLAGS) $(CPPFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

# -------------------- tar-ball --------------------
tar:
	@echo "creating zipped tar-ball ... "
	@cd .. ; tar -chvzf AsciiFromRoot/AsciiFromRoot.tar.gz AsciiFromRoot/AsciiFromRoot.cc AsciiFromRoot/RootUtilities.?? ;

# -------------------- clean --------------------

clean:
	rm  -f $(NAME) *.o
