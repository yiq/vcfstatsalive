#ifdef DUMMYSTATSCOLLECTOR_H
#define DUMMYSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"

namespace VcfStatsAlive {

	class DummyStatsCollector : public AbstractStatCollector {
		protected:
			virtual void processVariantImpl(const vcf::Variant& var) {;}
			virtual void appendJsonImpl(json_t * jsonRootObj) {;}
	};
}


#endif