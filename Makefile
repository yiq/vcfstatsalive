# Use custom to pass -DVCFLIB_PARITY to omit certain htslib API calls
CUSTOM=

CFLAGS=-g -std=c++11 $(CUSTOM)
INCLUDES=-Ilib/vcflib/src -Ilib/vcflib -Ilib/jansson-2.6/src -Ilib/htslib/

LDADDS=-lz -lstdc++ 
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin)
	LIBEX=dylib
else
	LDADDS+=-lhts
	LIBEX=so
endif

SOURCES=main.cpp \
		AbstractStatCollector.cpp \
		BasicStatsCollector.cpp 
PROGRAM=vcfstatsalive
PCH_SOURCE=vcfStatsAliveCommon.hpp
PCH=$(PCH_SOURCE).gch

PCH_FLAGS=-include $(PCH_SOURCE)

OBJECTS=$(SOURCES:.cpp=.o)

JANSSON=lib/jansson-2.6/src/.libs/libjansson.a
HTSLIB=lib/htslib/libhts.$(LIBEX)

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(PCH) $(OBJECTS) $(JANSSON) $(HTSLIB)
	$(CXX) $(CFLAGS)  $(HTSLIB) -v -o $@ $(OBJECTS) $(JANSSON) $(LDADDS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) $(PCH_FLAGS) -c $< -o $@

$(PCH): 
	$(CXX) $(CFLAGS) $(INCLUDES) -x c++-header $(PCH_SOURCE) -Winvalid-pch -o $@

.PHONY: $(PCH)

$(JANSSON):
	@if [ ! -d lib/jansson-2.6 ]; then cd lib; curl -o - http://www.digip.org/jansson/releases/jansson-2.6.tar.gz | tar -xzf - ; fi
	@cd lib/jansson-2.6; ./configure --disable-shared --enable-static; make; cd ../..

$(HTSLIB): 
	cd lib/htslib && autoheader 
	cd lib/htslib && autoconf
	cd lib/htslib && ./configure
	make -C lib/htslib libhts.$(LIBEX)
	make -C lib/htslib install-$(LIBEX)

clean:
	rm -rf $(OBJECTS) $(PROGRAM) $(PCH) *.dSYM

clean-dep:
	make -C lib/jansson-2.6 clean
	make -C lib/htslib clean


.PHONY: clean
