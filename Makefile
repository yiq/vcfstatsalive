INCLUDES=-Ivcflib/src -Ivcflib
LDADDS=-lz -lstdc++

SOURCES=main.cpp
STATIC_LIBS=vcflib/libvcf.a
PROGRAM=vcfstatsalive

DISORDER=vcflib/smithwaterman/disorder.c

OBJECTS=$(SOURCES:.cpp=.o)

all: $(PROGRAM)

.PHONY: all


$(PROGRAM): $(OBJECTS) $(STATIC_LIBS)
	$(CXX) -o $@ $(OBJECTS) $(STATIC_LIBS) $(DISORDER) $(LDADDS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

vcflib/libvcf.a:
	make -C vcflib libvcf.a

clean:
	rm -rf $(OBJECTS) $(PROGRAM)
	make -C vcflib clean

.PHONY: clean
