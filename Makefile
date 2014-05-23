INCLUDES=-Ivcflib
LDADDS=-lstdc++

SOURCES=main.cpp
STATIC_LIBS=vcflib/libvcf.a
PROGRAM=vcfstatsalive


OBJECTS=$(SOURCES:.cpp=.o)


all: $(PROGRAM)

.PHONY: all


$(PROGRAM): $(OBJECTS) $(STATIC_LIBS)
	$(CC) -o $@ $(OBJECTS) $(STATIC_LIBS) $(LDADDS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

vcflib/libvcf.a:
	make -C vcflib libvcf.a

clean:
	rm -rf $(OBJECTS) $(PROGRAM)
	make -C vcflib clean

.PHONY: clean
