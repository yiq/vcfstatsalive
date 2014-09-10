#include "BasicStatsCollector.h"

using namespace VcfStatsAlive;

static const unsigned int kCMTrailLength = 5;
static const double kCMThreshold = 0.001;

inline bool _isPurine(const string& allele) {
	if(allele.length() != 1) return false;
	if(allele == "A" || allele == "a" || allele == "G" || allele == "g") return true;
	return false;
}

inline bool _isPyrimidine(const string& allele) {
	if(allele.length() != 1) return false;
	if(allele == "C" || allele == "c" || allele == "T" || allele == "t") return true;
	return false;
}

BasicStatsCollector::BasicStatsCollector() :
	AbstractStatCollector(),
	_transitions(0),
	_transversions(0) {

	_stats.clear();

	_stats[kTotalRecords] = 0;
	_stats[kTsTvRatio] = 0;

#ifdef DEBUG
	StatMapT::iterator iter;
	for(iter = _stats.begin(); iter != _stats.end(); iter++) {
		std::cerr<<"Initializing: "<<iter->first<<std::endl;
	}
#endif

}

void BasicStatsCollector::processVariantImpl(const vcf::Variant& var) {
	// increment total variant counter
	++_stats[kTotalRecords];

	std::vector<std::string>::const_iterator altIter = var.alt.cbegin();
	for(;altIter != var.alt.cend();altIter++) {
		if(_isPurine(var.ref)) {
			if(_isPurine(*altIter)) _transitions++;
			else _transversions++;
		}
		else {
			if(_isPurine(*altIter)) _transversions++;
			else _transitions++;
		}
	}
}

void BasicStatsCollector::appendJsonImpl(json_t * jsonRootObj) {

	// update some stats
	_stats[kTsTvRatio] = double(_transitions) / double(_transversions);
	std::cerr<<_transitions<<", "<<_transversions<<std::endl;

	StatMapT::iterator sIter;
	for(sIter = _stats.begin(); sIter != _stats.end(); sIter++) {
		json_object_set_new(jsonRootObj, sIter->first.c_str(), json_real(sIter->second));
	}

	bool consensus = true;
}
