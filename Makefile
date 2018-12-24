# Use custom to pass -DVCFLIB_PARITY to omit certain htslib API calls
CUSTOM?=

CFLAGS=-g -std=c++11 $(CUSTOM)
INCLUDES=-Ilib/jansson-2.6/src -Ilib/htslib/

LDADDS=-lz -lm -lbz2 -llzma -lstdc++ -lcurl -lcrypto -lpthread
UNAME_S := $(shell uname -s)

SOURCES=main.cpp \
		AbstractStatCollector.cpp \
		BasicStatsCollector.cpp 
PROGRAM=vcfstatsalive vcfstats
PCH_SOURCE=vcfStatsAliveCommon.hpp
PCH=$(PCH_SOURCE).gch

PCH_FLAGS=-include $(PCH_SOURCE)

OBJECTS=AbstractStatCollector.o \
	BasicStatsCollector.o 

JANSSON=lib/jansson-2.6/src/.libs/libjansson.a
HTSLIB=$(HTSLIB_HOME)/lib/libhts.a

all: $(PROGRAM)

.PHONY: all

vcfstatsalive: $(PCH) $(OBJECTS) main.cpp $(JANSSON)
	$(CXX) $(CFLAGS) $(PCH_FLAGS) -v -o $@ $(OBJECTS) main.cpp $(JANSSON) $(HTSLIB) $(LDADDS)

vcfstats: $(PCH) $(OBJECTS) vcfstats.cpp $(JANSSON)
	$(CXX) $(CFLAGS) $(PCH_FLAGS) -v -o $@ $(OBJECTS) vcfstats.cpp $(JANSSON) $(HTSLIB) $(LDADDS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) $(PCH_FLAGS) -c $< -o $@

$(PCH): 
	$(CXX) $(CFLAGS) $(INCLUDES) -x c++-header $(PCH_SOURCE) -Winvalid-pch -o $@

.PHONY: $(PCH)

$(JANSSON):
	@if [ ! -d lib/jansson-2.6 ]; then cd lib; curl -o - http://www.digip.org/jansson/releases/jansson-2.6.tar.gz | tar -xzf - ; fi
	@cd lib/jansson-2.6; ./configure --disable-shared --enable-static; make; cd ../..

clean:
	rm -rf $(OBJECTS) $(PROGRAM) $(PCH) *.dSYM

clean-dep:
	make -C lib/jansson-2.6 clean
	make -C lib/htslib clean


.PHONY: clean
