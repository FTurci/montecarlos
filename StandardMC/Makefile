##########################
# Macros: name = data
##########################

CXX=g++-6
CXXFLAGS= -std=c++11 -Wall -Wno-deprecated -O3 
CCFLAGS = $(CXXFLAGS)  -lgsl -lgslcblas
LIBS= -lm
LDFLAGS =

#########################
#Rules: target: source
#		command
#########################

#.SUFFIXES: .cc .c .C

#.C.o:
#	$(CC) 	-c 	$(CFLAGS) 	-o 	$@ 	$<

all: 	monte

monte:	 liquid.o main.o compass_search.o
	$(CXX) $(LDFLAGS) $(LIBS) liquid.o  main.o  compass_search.o -o monte

clean: 
	rm  liquid.o  main.o compass_search.o
