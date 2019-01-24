#include "BasicStatsCollector.h"

#include <cmath>

using namespace std;

namespace VcfStatsAlive {
	inline unsigned int StringToUInt(const string &text)
	{
		std::stringstream ss(text);
		unsigned int result;
		return ss >> result ? result : 0;
	}

	inline double StringToDouble(const string &text)
	{
		std::stringstream ss(text);
		double result;
		return ss >> result ? result : 0;
	}

	template <typename T>
	inline int valueToBinCloseInterval(T value, T minValue, T maxValue, int bins) {
		return  floor((value - minValue) / (maxValue - minValue) * bins);
	}

	template <typename T>
	inline int valueToBinOpenInterval(T value, T minValue, T maxValue, int bins) {
		if (value < minValue) return 0;
		else return valueToBinCloseInterval(value, minValue, maxValue, bins-1) + 1;
	}
}

using namespace VcfStatsAlive;

static const unsigned int kCMTrailLength = 5;
static const double kCMThreshold = 0.001;
static const double kLogAFLowerBound = -5.0;
static const double kLogAFUpperBound = 0.0;

inline bool _isPurine(const char allele) {
	//if(allele.length() != 1) return false;
	if(allele == 'A' || allele == 'a' || allele == 'G' || allele == 'g') return true;
	return false;
}

inline bool _isPyrimidine(const char allele) {
	//if(allele.length() != 1) return false;
	if(allele == 'C' || allele == 'c' || allele == 'T' || allele == 't') return true;
	return false;
}

inline size_t base2Idx(const char base) {
	if(base == 'A' || base == 'a') return 0;
	else if (base == 'G' || base == 'g') return 1;
	else if (base == 'C' || base == 'c') return 2;
	else if (base == 'T' || base == 't') return 3;
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

BasicStatsCollector::BasicStatsCollector(int qualLower, int qualUpper, bool logScaleAF) :
	AbstractStatCollector(),
	_transitions(0),
	_transversions(0),
	usingLogScaleAF(logScaleAF),
	kQualHistLowerbound(qualLower),
	kQualHistUpperbound(qualUpper),
	m_qualityDist(qualUpper - qualLower + 1 + 2, 0) {

	_stats.clear();

	_stats[kTotalRecords] = 0;
	_stats[kTsTvRatio] = 0;

	_alleleFreqBins = logScaleAF ? 52 : 51;
	m_alleleFreqHist = (unsigned int *)malloc(sizeof(unsigned int) * _alleleFreqBins);


	memset(m_alleleFreqHist, 0, sizeof(unsigned int) * _alleleFreqBins);
	memset(m_mutationSpec, 0, sizeof(unsigned int) * 4 * 4);
	memset(m_variantTypeDist, 0, sizeof(unsigned int) * static_cast<unsigned int>(VT_SIZE));


#ifdef DEBUG
	StatMapT::iterator iter;
	for(iter = _stats.begin(); iter != _stats.end(); iter++) {
		std::cerr<<"Initializing: "<<iter->first<<std::endl;
	}
#endif

}

BasicStatsCollector::~BasicStatsCollector() {
	free(m_alleleFreqHist);
}

void BasicStatsCollector::updateTsTvRatio(bcf1_t* var, int altIndex, bool isSnp) {
	// TsTv Ratio - Only evaluate SNPs 

	char ref = var->d.allele[0][0];
	char htsAlt = var->d.allele[altIndex][0];
#ifdef VCFLIB_PARITY
	if(var->d.allele[0][1] == 0 && var->d.allele[altIndex][1] == 0 && var->d.allele[0][0] != '.' && var->d.allele[altIndex][0] != '.')
#else
	if(isSnp) 
#endif
    {
		if(_isPurine(ref)) {
			if(_isPurine(htsAlt)) {
				_transitions++;
			}
			else if(_isPyrimidine(htsAlt)){
				_transversions++;
			} 
		}
		else {
			if(_isPurine(htsAlt)) {
				_transversions++;
			}
			else if(_isPyrimidine(htsAlt)){
				_transitions++;
			} 
		}	
	}

}

void BasicStatsCollector::updateMutationSpectrum(bcf1_t* var, int altIndex, bool isSnp) {
	// Mutation Spectrum
#ifdef VCFLIB_PARITY
	if(var->d.allele[0][1] != 0 || var->d.allele[altIndex][1] != 0) return;
#else
	if(!isSnp) return;
#endif

	size_t firstIdx = base2Idx(var->d.allele[0][0]);
	size_t secondIdx = base2Idx(var->d.allele[altIndex][0]);

	if(firstIdx < 4 && secondIdx < 4) {
		m_mutationSpec[firstIdx][secondIdx]++;
	}
}

void BasicStatsCollector::updateAlleleFreqHist(bcf_hdr_t* hdr, bcf1_t* var) {
	// Allele Frequency Histogram
	int alleleFreqBin;
	double alleleFreq = 0;
	int count = 0;

	float* pfAlleleFreq = NULL;
	if(bcf_get_info_float(hdr, var, "AF", &pfAlleleFreq, &count) >= 0 && count > 0) {
		alleleFreq = *pfAlleleFreq;
	}
	else {
		int* piDepth = NULL;
		int* piRefObsrv = NULL;
		count = 0;
		if(bcf_get_info_int32(hdr, var, "DP", &piDepth, &count) >= 0 && count > 0 && 
		   bcf_get_info_int32(hdr, var, "RO", &piRefObsrv, &count) >= 0 && count > 0) {
					alleleFreq = ( *piDepth - *piRefObsrv ) / ((double)*piDepth);
		}

		free(piRefObsrv);
		free(piDepth);
	}

	free(pfAlleleFreq);

	if(alleleFreq == 0) return;

	if (usingLogScaleAF) {
		double logAF = log10(alleleFreq);
		alleleFreqBin = (int) valueToBinOpenInterval(logAF, kLogAFLowerBound, kLogAFUpperBound, _alleleFreqBins - 1);
	}
	else {
		alleleFreqBin = (int) valueToBinCloseInterval(alleleFreq, 0.0, 1.0, _alleleFreqBins - 1);
	}

	//assert(alleleFreq <= 1);
	//assert(alleleFreq >= 0);
	assert(alleleFreqBin >= 0);
	assert(alleleFreqBin < _alleleFreqBins);

	m_alleleFreqHist[alleleFreqBin]++;
}

void BasicStatsCollector::updateVariantTypeDist(bcf1_t* var, int altIndex, int refLength) {
	// Type Distribution
	VariantTypeT vt;
	int altLength = strlen(var->d.allele[altIndex]);

	// Do not use bcf_is_snp here because it enforces its logic across all alternates.
	if(refLength == 1 && altLength == 1) {
		vt = VT_SNP;
	}
    else if (refLength == altLength) {
        vt = VT_MNP;
    }
    else if (strstr(var->d.allele[altIndex], "<") != nullptr) {
        vt = VT_SV;
    }
	else if (altLength > refLength) {
		vt = VT_INS;
		updateIndelSizeDist(refLength, altLength);
	}
	else if (altLength < refLength) {
		vt = VT_DEL;
		updateIndelSizeDist(refLength, altLength);
	}
	else {
		vt = VT_OTHER;
	}

	m_variantTypeDist[static_cast<unsigned int>(vt)]++;
}

void BasicStatsCollector::updateQualityDist(float qual) {

	size_t bin = 0;

	int intQual = int(qual);

	if (intQual < kQualHistLowerbound)
		bin = (kQualHistUpperbound - kQualHistLowerbound + 1);
	else if (intQual > kQualHistUpperbound + 1)
		bin = (kQualHistUpperbound - kQualHistLowerbound + 2);
	else
		bin = (intQual - kQualHistLowerbound);

	m_qualityDist[bin] += 1;

}

void BasicStatsCollector::updateIndelSizeDist(int refLength, int altLength) {

	long indelSize = long(altLength) - long(refLength);
	if(m_indelSizeDist.find(indelSize) == m_indelSizeDist.end())
		m_indelSizeDist[indelSize] = 1;
	else
		m_indelSizeDist[indelSize] += 1;
}

void BasicStatsCollector::processVariantImpl(bcf_hdr_t* hdr, bcf1_t* var) {
	// increment total variant counter
	++_stats[kTotalRecords];

	bool isSnp = bcf_is_snp(var);
	int refLength = strlen(var->d.allele[0]);

	for(int altIndex = 1; altIndex < var->n_allele; altIndex++) {
		updateTsTvRatio(var, altIndex, isSnp);
		updateMutationSpectrum(var, altIndex, isSnp);
		updateVariantTypeDist(var, altIndex, refLength);
	}

	updateAlleleFreqHist(hdr, var);
	updateQualityDist(var->qual);

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
	json_t * j_af_hist_bins = json_object();
	for(size_t i=0; i<50; i++) {
		if (m_alleleFreqHist[i] > 0) {
			std::stringstream labelSS; labelSS << i;
			json_object_set_new(j_af_hist_bins, labelSS.str().c_str(), json_integer(m_alleleFreqHist[i]));
		}
	}

	json_object_set_new(j_af_hist, "usingLogScaleAF", json_boolean(usingLogScaleAF));
	if (usingLogScaleAF) {
		json_object_set_new(j_af_hist, "logAFHistLowerBound", json_real(kLogAFLowerBound));
		json_object_set_new(j_af_hist, "logAFHistUpperBound", json_real(kLogAFUpperBound));
	}
	json_object_set_new(j_af_hist, "afHistBins", j_af_hist_bins);
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
			case VT_MNP:
				label = "MNP";
                break;
			case VT_SV:
				label = "SV";
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

	// Quality Distribution
	json_t * j_qual_dist = json_object();
	json_object_set_new(j_qual_dist, "qualHistLowerBound", json_integer(kQualHistLowerbound));
	json_object_set_new(j_qual_dist, "qualHistUpperBound", json_integer(kQualHistUpperbound));
	json_t *j_qual_dist_bin = json_object();
	for(size_t i=0; i<m_qualityDist.size() - 2; i++) {
		if (m_qualityDist[i] == 0) continue;
		std::stringstream labelSS; labelSS << i;
		json_object_set_new(j_qual_dist_bin, labelSS.str().c_str(), json_integer(m_qualityDist[i]));
	}
	json_object_set_new(j_qual_dist, "regularBins", j_qual_dist_bin);
	json_object_set_new(j_qual_dist, "lowerBin", json_integer(m_qualityDist[kQualHistUpperbound - kQualHistLowerbound + 1]));
	json_object_set_new(j_qual_dist, "upperBin", json_integer(m_qualityDist[kQualHistUpperbound - kQualHistLowerbound + 2]));

	json_object_set_new(jsonRootObj, "qual_dist", j_qual_dist);

	// Indel Size Dist
	json_t * j_indel_size= json_object();
	for(map<long, size_t>::iterator it = m_indelSizeDist.begin(); it != m_indelSizeDist.end(); it++) {
		std::stringstream labelSS; labelSS << it->first;
		json_object_set_new(j_indel_size, labelSS.str().c_str(), json_integer(it->second));
	}
	json_object_set_new(jsonRootObj, "indel_size", j_indel_size);



	bool consensus = true;
}
