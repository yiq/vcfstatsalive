#include "BasicStatsCollector.h"

#include <cmath>

namespace VcfStatsAlive {
	unsigned int StringToUInt(const string &text)
	{
		std::stringstream ss(text);
		unsigned int result;
		return ss >> result ? result : 0;
	}

	double StringToDouble(const string &text)
	{
		std::stringstream ss(text);
		double result;
		return ss >> result ? result : 0;
	}
}

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

inline size_t base2Idx(const string& base) {
	if(base == "A" || base == "a") return 0;
	else if (base == "G" || base == "g") return 1;
	else if (base == "C" || base == "c") return 2;
	else if (base == "T" || base == "t") return 3;
	else return 4;
}

inline char idx2Base(size_t idx) {
	switch(idx) {
		case 0:
			return 'A';
		case 1:
			return 'G';
		case 2:
			return 'C';
		case 3:
			return 'T';
		default:
			return ' ';
	}
}

BasicStatsCollector::BasicStatsCollector() :
	AbstractStatCollector(),
	_transitions(0),
	_transversions(0) {

	_stats.clear();

	_stats[kTotalRecords] = 0;
	_stats[kTsTvRatio] = 0;

	memset(m_alleleFreqHist, 0, sizeof(unsigned int) * 50);
	memset(m_mutationSpec, 0, sizeof(unsigned int) * 4 * 4);
	memset(m_variantTypeDist, 0, sizeof(unsigned int) * VT_SIZE);

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

	std::vector<std::string>::const_iterator altIter = var.alt.begin();
	for(;altIter != var.alt.end();altIter++) {
		// TsTv Ratio
		if(_isPurine(var.ref)) {
			if(_isPurine(*altIter)) {
				_transitions++;
			}
			else {
				_transversions++;
			}
		}
		else {
			if(_isPurine(*altIter)) {
				_transversions++;
			}
			else {
				_transitions++;
			}
		}

		// Mutation Spectrum
		size_t firstIdx = base2Idx(var.ref);
		size_t secondIdx = base2Idx(*altIter);

		if(firstIdx < 4 && secondIdx < 4) {
			m_mutationSpec[firstIdx][secondIdx]++;
		}

		// Type Distribution
		int vt=-1;
		if(var.ref.size() == 1 && altIter->size() == 1) {
			vt = VT_SNP;
		}
		else if (var.ref.size() == 1 && altIter->size() > 1) {
			vt = VT_INS;
		}
		else if (var.ref.size() > 1 && altIter->size() == 1) {
			vt = VT_DEL;
		}
		else {
			vt = VT_OTHER;
		}

		m_variantTypeDist[vt]++;
	}

	// Allele Frequency Histogram
	int alleleFreqBin;
	double alleleFreq;

	if(var.info.find("AF") != var.info.end()) {
		alleleFreq = StringToDouble(var.info.at("AF")[0]);
	}
	else {
		unsigned int depth = StringToUInt(var.info.at("DP")[0]);
		unsigned int refObsrv = StringToUInt(var.info.at("RO")[0]);
		double alleleFreq = ( depth - refObsrv ) / ((double)depth);
	}

	alleleFreqBin = (int) ceil(alleleFreq/2.0*100.0) - 1;

	if(alleleFreqBin == -1) alleleFreqBin = 0;

	assert(alleleFreq <= 1);
	assert(alleleFreq >= 0);
	
	m_alleleFreqHist[alleleFreqBin]++;
}

void BasicStatsCollector::appendJsonImpl(json_t * jsonRootObj) {

	// update some stats
	_stats[kTsTvRatio] = double(_transitions) / double(_transversions);

	StatMapT::iterator sIter;
	for(sIter = _stats.begin(); sIter != _stats.end(); sIter++) {
		json_object_set_new(jsonRootObj, sIter->first.c_str(), json_real(sIter->second));
	}

	// Allele Frequency Histogram
	json_t * j_af_hist = json_object();
	for(size_t i=0; i<50; i++) {
		if (m_alleleFreqHist[i] > 0) {
			std::stringstream labelSS; labelSS << i;
			json_object_set_new(j_af_hist, labelSS.str().c_str(), json_integer(m_alleleFreqHist[i]));
		}
	}
	json_object_set_new(jsonRootObj, "af_hist", j_af_hist);

	// Mutation spectrum
	json_t * j_mut_spec = json_object();
	for(size_t first=0; first<4; first++) {
		std::stringstream labelSS; labelSS << idx2Base(first);
		json_t * j_spec_array = json_array();
		for(size_t second = 0; second<4; second++) {
			json_array_append_new(j_spec_array, json_integer(m_mutationSpec[first][second]));
		}
		json_object_set_new(j_mut_spec, labelSS.str().c_str(), j_spec_array);
	}
	json_object_set_new(jsonRootObj, "mut_spec", j_mut_spec);

	// Mutation type
	json_t * j_mut_type = json_object();
	for(size_t vt = 0; vt < VT_SIZE; vt++) {
		string label;
		switch(vt) {
			case VT_SNP:
				label = "SNP";
				break;
			case VT_INS:
				label = "INS";
				break;
			case VT_DEL:
				label = "DEL";
				break;
			default:
				label = "OTHER";
				break;
		}

		json_object_set_new(j_mut_type, label.c_str(), json_integer(m_variantTypeDist[vt]));
	}
	json_object_set_new(jsonRootObj, "var_type", j_mut_type);

	bool consensus = true;
}
