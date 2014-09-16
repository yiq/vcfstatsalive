#ifndef BASICSTATSCOLLECTOR_H
#define BASICSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"

namespace VcfStatsAlive {

	static std::string const kTotalRecords = "TotalRecords";
	static std::string const kTsTvRatio = "TsTvRatio";

	typedef std::map<std::string, double> StatMapT;

	class BasicStatsCollector : public AbstractStatCollector {

		protected:
			StatMapT _stats;

			size_t _transitions;
			size_t _transversions;
			unsigned int m_alleleFreqHist[50];
			unsigned int m_mutationSpec[4][4];

			virtual void processVariantImpl(const vcf::Variant& var);
			virtual void appendJsonImpl(json_t * jsonRootObj);

		public:
			BasicStatsCollector();
	};
}

#endif
