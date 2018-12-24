/* ByGenotypeStratifier.h
 *
 * Collect stats while stratifing over genotype
 *
 * Because the stratifier will need to be able to create
 * new stats collectors of the correct type, this class
 * is defined as a template. User should instanciate this class
 * using the desired stats collector type
 */

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

                auto gt_cat = genotype_category(gt_arr, ngt_arr);

                if(gt_cat == "REF") return;

                auto coll = m_collectors.find(gt_cat);

                if(coll == m_collectors.end()) {
                    std::cerr<<"collector created for "<<gt_cat<<std::endl;
                    m_collectors[gt_cat] = new CollectorT();
                }

                m_collectors[gt_cat]->processVariant(hdr, var);
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

            // convert htslib genotype array to category string
            std::string genotype_category(const int32_t *gt_arr, const int32_t ngt_arr) {
                /* categories:
                 *
                 * 1. heterozygous.     E.g. 0/1, 1/2, 2|1
                 * 2. homozygous alt.   E.g. 1/1, 2|2
                 */

                auto gt1 = (gt_arr[0] >> 1) - 1;
                auto gt2 = (gt_arr[1] >> 1) - 1;

                if(gt1 < 0 || gt2 < 0) return "MISSING";
                if(gt1 == 0 && gt2 == 0) return "REF";
                if(gt1 == gt2) return "HOMO";
                return "HET";
           }
    };
}

#endif
