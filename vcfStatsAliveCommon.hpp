/**
 * @file bamstatsAliveCommon.hpp
 * This is the common header for the entire project.
 *
 * This header file will be compiled into precompiled header, and
 * automatically be included in all source code's compilation.
 *
 * @author Yi Qiao
 */

#ifndef VCFSTATSALIVECOMMON_H
#define VCFSTATSALIVECOMMON_H

// Standard C/C++ Headers
#include <iostream>
#include <sstream>
#include <cstdlib>

// Commonly used STL classes
#include <string>
#include <vector>
#include <map>
#include <algorithm>

// Utilize assertions
#include <assert.h>

// Jansson JSON manipulation library
#include <jansson.h>

// Include htslib
namespace htslib
{
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
}

#define kstring_t htslib::kstring_t
#include "Variant.h"


// Logging facility for Debug
#ifdef RELEASE
namespace VcfStatsAlive {
	namespace Log {
		static std::ostream bitBucket(0);
	}
}
#undef LOGS
#define LOGS (Log::bitBucket)
#endif

#ifdef DEBUG
#undef LOGS
#define LOGS (std::cerr)
#endif

#endif
