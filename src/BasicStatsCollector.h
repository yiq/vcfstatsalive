#ifndef BASICSTATSCOLLECTOR_H
#define BASICSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"

namespace VcfStatsAlive {

	static std::string const kTotalRecords = "TotalRecords";
	static std::string const kTsTvRatio = "TsTvRatio";

	using StatMapT = std::map<std::string, double>;

	typedef enum {
		VT_SNP = 0,
		VT_INS,
		VT_DEL,
        VT_MNP,
        VT_SV,
		VT_OTHER,
		VT_SIZE
	} VariantTypeT;

	class BasicStatsCollector : public AbstractStatCollector {

		protected:
			StatMapT _stats;

			size_t _transitions;
			size_t _transversions;

			const int kQualHistLowerbound;
			const int kQualHistUpperbound;

			unsigned int *m_alleleFreqHist;
			size_t _alleleFreqBins;
			bool usingLogScaleAF;
			std::vector<int> m_qualityDist;
			unsigned int m_mutationSpec[4][4];
			unsigned int m_variantTypeDist[VT_SIZE];
			std::map<long, size_t> m_indelSizeDist;


			virtual void processVariantImpl(bcf_hdr_t* hdr, bcf1_t* var) override;
			virtual void appendJsonImpl(json_t * jsonRootObj) override;

            void updateTsTvRatio(bcf1_t* var, int altIndex, bool isSnp);
            void updateMutationSpectrum(bcf1_t* var, int altIndex, bool isSnp);
            void updateAlleleFreqHist(bcf_hdr_t* hdr, bcf1_t* var);
            void updateQualityDist(float qual);
            void updateVariantTypeDist(bcf1_t* var, int altIndex, int refLength);
            void updateIndelSizeDist(int refLength, int altLength);

		public:
			BasicStatsCollector(int qualLower, int qualUpper, bool logScaleAF = false);
			virtual ~BasicStatsCollector();
	};
}

#endif
