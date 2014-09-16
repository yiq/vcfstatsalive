CFLAGS=-g
INCLUDES=-Ilib/vcflib/src -Ilib/vcflib -Ilib/jansson-2.6/src
LDADDS=-lz -lstdc++

SOURCES=main.cpp \
		AbstractStatCollector.cpp \
		BasicStatsCollector.cpp 
PROGRAM=vcfstatsalive
PCH_SOURCE=vcfStatsAliveCommon.hpp
PCH=$(PCH_SOURCE).gch

PCH_FLAGS=-include $(PCH_SOURCE)

OBJECTS=$(SOURCES:.cpp=.o)

JANSSON=lib/jansson-2.6/src/.libs/libjansson.a
VCFLIB=lib/vcflib/libvcf.a
DISORDER=lib/vcflib/smithwaterman/disorder.c

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(PCH) $(OBJECTS) $(VCFLIB) $(JANSSON)
	$(CXX) $(CFLAGS) -o $@ $(OBJECTS) $(VCFLIB) $(JANSSON) $(DISORDER) $(LDADDS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) $(PCH_FLAGS) -c $< -o $@

$(PCH): 
	$(CXX) $(CFLAGS) $(INCLUDES) -x c++-header $(PCH_SOURCE) -Winvalid-pch -o $@

.PHONY: $(PCH)

$(JANSSON):
	@if [ ! -d lib/jansson-2.6 ]; then cd lib; curl -o - http://www.digip.org/jansson/releases/jansson-2.6.tar.gz | tar -xzf - ; fi
	@cd lib/jansson-2.6; ./configure --disable-shared --enable-static; make; cd ../..

$(VCFLIB):
	make -C lib/vcflib libvcf.a

clean:
	rm -rf $(OBJECTS) $(PROGRAM) $(PCH) *.dSYM

clean-dep:
	make -C lib/vcflib clean
	make -C lib/jansson-2.6 clean

.PHONY: clean
