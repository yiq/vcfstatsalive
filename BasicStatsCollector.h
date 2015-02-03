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

			unsigned int m_alleleFreqHist[50];
			vector<int> m_qualityDist;
			unsigned int m_mutationSpec[4][4];
			unsigned int m_variantTypeDist[VT_SIZE];
			map<long, size_t> m_indelSizeDist;


			virtual void processVariantImpl(const vcf::Variant& var) override;
			virtual void appendJsonImpl(json_t * jsonRootObj) override;

		public:
			BasicStatsCollector(int qualLower, int qualUpper);

		private:
			void updateTsTvRatio(const vcf::Variant& var, const string& alt);
			void updateMutationSpectrum(const vcf::Variant& var, const string& alt);
			void updateAlleleFreqHist(const vcf::Variant& var);
			void updateQualityDist(const vcf::Variant& var);
			void updateVariantTypeDist(const vcf::Variant& var, const string& alt);
			void updateIndelSizeDist(const vcf::Variant& var, const string& alt);
	};
}

#endif
