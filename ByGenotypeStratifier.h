#ifndef BYGENOTYPESTRATIFIER_H
#define BYGENOTYPESTRATIFIER_H

#pragma once

#include "AbstractStatCollector.h"
#include <sstream>

namespace VcfStatsAlive {

    template <class CollectorT>
    class ByGenotypeStratifier : public AbstractStatCollector {
        protected:
            virtual void processVariantImpl(bcf_hdr_t* hdr, bcf1_t* var) override {
                int32_t ngt, *gt_arr=nullptr, ngt_arr=0;
                ngt = bcf_get_genotypes(hdr, var, &gt_arr, &ngt_arr);
                std::unique_ptr<int32_t> uniq_gt_arr(gt_arr);

                auto gt_str = genotype_str(gt_arr, ngt_arr);

                auto coll = m_collectors.find(gt_str);

                if(coll == m_collectors.end()) {
                    std::cout<<"collector created for "<<gt_str<<std::endl;
                    m_collectors[gt_str] = new CollectorT();
                }

                m_collectors[gt_str]->processVariant(hdr, var);
            }
            
            virtual void appendJsonImpl(json_t* jsonRootObj) override {
                for(auto& gt : m_collectors) {
                    json_t * j_gt = json_object();
                    gt.second->appendJson(j_gt);
                    json_object_set_new(jsonRootObj, gt.first.c_str(), j_gt);
                }
            }

        public:
            ByGenotypeStratifier() : AbstractStatCollector() { }
            virtual ~ByGenotypeStratifier() {
                for(auto& gt : m_collectors) delete gt.second;
            }

        private:
            std::map<std::string, CollectorT*> m_collectors;

            // convert htslib genotype array to vcf genotype string
            std::string genotype_str(const int32_t *gt_arr, const int32_t ngt_arr) {
                std::stringstream ss;
                for(int32_t i=0; i<ngt_arr; i++) {
                    if(i > 0) ss<<(gt_arr[i] & 1 ? '|':'/');
                    auto gt = (gt_arr[i]>>1);
                    ss<< (gt > 0 ? (gt - 1) : '.');
                }

                return ss.str();
            }
    };
}

#endif
