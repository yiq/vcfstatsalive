INCLUDES=-Ivcflib


SOURCES=main.cpp
STATIC_LIBS=vcflib/libvcf.a
PROGRAM=vcfstatsalive


OBJECTS=$(SOURCES:.cpp=.o)


all: $(PROGRAM)

.PHONY: all


$(PROGRAM): $(OBJECTS) $(STATIC_LIBS)
	$(CC) $? $(STATIC_LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

libvcf/libvcf.a:
	make -C vcflib libvcf.a

clean:
	rm -rf $(OBJECTS) $(PROGRAM)

.PHONY: clean