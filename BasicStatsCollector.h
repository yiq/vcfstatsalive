#ifndef BASICSTATSCOLLECTOR_H
#define BASICSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"
#include "AbstractChangeMonitor.h"

namespace VcfStatsAlive {

	static std::string const kTotalRecords = "TotalRecords";
	static std::string const kTsTvRatio = "TsTvRatio";

	typedef std::map<std::string, unsigned int> StatMapT;
	typedef std::map<std::string, AbstractChangeMonitor<double> *> ChangeMonitorMapT;

	class BasicStatsCollector : public AbstractStatCollector {

		protected:
			StatMapT _stats;
			ChangeMonitorMapT _monitors;

			size_t _transitions;
			size_t _transversions;

			virtual void processVariantImpl(const vcf::Variant& var);
			virtual void appendJsonImpl(json_t * jsonRootObj);

		public:
			BasicStatsCollector();
	};
}

#endif
