CC=clang
CFLAGS=-c -Wall -O2 -std=c++0x
LDFLAGS=
CPP_FILES=$(wildcard src/*.cc) 
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cc=.o)))
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=align

$(EXECUTABLE): $(OBJ_FILES)
	g++ $(LDFLAGS) -o $@ $^

obj/%.o: src/%.cc
	g++ $(CFLAGS) -c -o $@ $<

clean:
	rm obj/*.o 
