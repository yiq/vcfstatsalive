#include "AbstractStatCollector.h"

using namespace VcfStatsAlive;

AbstractStatCollector::AbstractStatCollector(const std::string* sampleName) {
	_children.clear();
}

AbstractStatCollector::~AbstractStatCollector() {

}

void AbstractStatCollector::addChild(StatCollectorPtr child) {

	// Make sure that the input is good
	if(child.get() == nullptr) return;

	// Make sure the input is not yet a child
	auto loc = std::find(_children.cbegin(), _children.cend(), child);
	if(loc != _children.cend()) return;

	// Insert the input to the end of the children list
	_children.push_back(child);
}

void AbstractStatCollector::removeChild(StatCollectorPtr child) {

	// Make sure that the input is good
	if(child.get() == nullptr) return;

	// Make sure the input is a child
	auto loc = std::find(_children.begin(), _children.end(), child);
	if(loc == _children.end()) return;

	_children.erase(loc);
}

void AbstractStatCollector::processVariant(const vcf::Variant& var, htslib::bcf_hdr_t* hdr, htslib::bcf1_t* htsVar) {

	this->processVariantImpl(var, hdr, htsVar);

	for(auto iter = _children.begin(); iter != _children.end(); iter++) {
		(*iter)->processVariant(var, hdr, htsVar);
	}
}

json_t * AbstractStatCollector::appendJson(json_t * jsonRootObj) {
	if(jsonRootObj == NULL)
		return NULL;

	this->appendJsonImpl(jsonRootObj);
	
	for(auto iter = _children.begin(); iter != _children.end(); iter++) {
		(*iter)->appendJson(jsonRootObj);
	}

	return jsonRootObj;
}

bool AbstractStatCollector::isSatisfiedImpl() {
	return false;
}

bool AbstractStatCollector::isSatisfied() {
	if (not this->isSatisfiedImpl()) return false;

	bool isChildrenSatisfied = true;

	for(auto iter = _children.begin(); iter != _children.end(); iter++) {
		isChildrenSatisfied = isChildrenSatisfied && (*iter)->isSatisfied();
	}

	return isChildrenSatisfied;
}
