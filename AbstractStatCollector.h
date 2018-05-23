#ifndef ABSTRACTSTATCOLLECTOR_H
#define ABSTRACTSTATCOLLECTOR_H

#pragma once

#include <memory>

namespace VcfStatsAlive {

	class AbstractStatCollector;
	
	using StatCollectorPtr = std::shared_ptr<AbstractStatCollector>;

	using StatCollectorPtrVec = std::vector<StatCollectorPtr>;

	/**
	 * The base class for all statistics collectors
	 *
	 * A statistics collector will implement two virtual functions: 
	 *   - processVariant() to update statistics
	 *   - appendJson() to create the json representation of the statistics
	 *
	 * These statistics collectors can be organized into a tree with the
	 * addChild() and removeChild() functions. User code will only need to call
	 * the public processVariant() and appendJson() functions on the root
	 * object, and the action will be propagated across all child nodes. The
	 * actual implementation of specific collectors is encapsulated by the
	 * protected processVariantImpl() and appendJsonImpl() functions
	 */
	class AbstractStatCollector {
		protected:
			StatCollectorPtrVec _children;

			/**
			 * Process the variant and update statistics
			 *
			 * @param hdr The vcf file header information
			 * @param var The hstlib variant record
			 */
			virtual void processVariantImpl(bcf_hdr_t* hdr, bcf1_t* var) = 0;

			/**
			 * Append statistics as json
			 *
			 * @param jsonRootObj The json root object to which the outputs are appended
			 */
			virtual void appendJsonImpl(json_t * jsonRootObj) = 0;

			/** 
			 * Check if the statistics collector is satisfied with the data it
			 * has seen so far. Note that the defualt implementation of this
			 * function always returns false, which means that by default, a
			 * collector will keep processing reads indefinitely.
			 *
			 * @return true if the data is considered to be sufficient, false otherwise
			 */
			virtual bool isSatisfiedImpl();

		public:
			AbstractStatCollector(const std::string* sampleName = NULL);
			virtual ~AbstractStatCollector();

			/**
			 * Add a statistics collector as the child of the current collector
			 * 
			 * @param child The child collector to be added
			 */
			void addChild(StatCollectorPtr child);

			/**
			 * Remove a statistics collector from the children list of the current collector
			 *
			 * @param child The child collector to be removed
			 */
			void removeChild(StatCollectorPtr child);

			/**
			 * Process an variant by the collector tree
			 *
			 * The variant will be passed to the processVariantImpl function of the
			 * current collector, and the processVariant function of all the children
			 * collectors.
			 *
			 * @param hdr The vcf file header information
			 * @param var The htslib variant 
			 */
			void processVariant(bcf_hdr_t* hdr, bcf1_t* var);

			/**
			 * Create json of the collector tree
			 *
			 * Json objects representing all the statistics collected by the
			 * current tree will be appended into the given json root object.
			 * If no root object is given, a new one will be created. Current
			 * collector's appendJsonImpl function and all the children's 
			 * appendJson function will be called on the root object.
			 *
			 * @param jsonRootObj The json root object to which all statistics are appended to
			 * @return The json object representing the root of all the statistics
			 */
			json_t * appendJson(json_t * jsonRootObj = NULL);

			/**
			 * Check satisfy-ness of the collector tree
			 *
			 * @return true if all collectors in the tree are satisfied, false otherwise
			 */
			bool isSatisfied();
	};

}

#endif
