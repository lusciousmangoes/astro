PROGRAM = dnc
OBJECTS = dnc.o build_tree.o calc_forces.o particle_io.o move.o

CFLAGS = -Wall -O4
LDLIBS = -lm

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)

install:
	cp dnc ../bin/dnc

clean:
	rm -f $(PROGRAM) $(OBJECTS)
