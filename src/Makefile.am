if DEBUG
AM_CFLAGS = -g -O0
AM_CXXFLAGS = -g -O0
else
AM_CFLAGS = -O2
AM_CXXFLAGS = -O2
endif

if PROFILING
AM_CFLAGS += -pg
AM_CXXFLAGS += -pg
endif

AM_CPPFLAGS = -include vcfStatsAliveCommon.hpp

bin_PROGRAMS=vcfstatsalive vcfstats

vcfstatsalive_SOURCES = main.cpp AbstractStatCollector.cpp BasicStatsCollector.cpp
vcfstats_SOURCES = vcfstats.cpp AbstractStatCollector.cpp BasicStatsCollector.cpp
