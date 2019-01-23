/* BySampleStratifier.h
 *
 * Collect stats while stratifing over samples
 *
 * Because the stratifier will need to be able to create
 * new stats collectors of the correct type, this class
 * is defined as a template. User should instanciate this class
 * using the desired stats collector type
 */

#ifndef BYSAMPLESTRATIFIER_H
#define BYSAMPLESTRATIFIER_H

#pragma once

#include "AbstractStatCollector.h"

namespace VcfStatsAlive {

    template <class CollectorT>
    class BySampleStratifier : public AbstractStatCollector {
        protected:
            virtual void processVariantImpl(bcf_hdr_t* hdr, bcf1_t *var) override {

                bcf1_t* subset_rec; // temporary bcf record to hold subsetted record

                for(int sample_idx = 0; sample_idx < hdr->n[BCF_DT_SAMPLE]; sample_idx++) {
                    auto* sample = hdr->id[BCF_DT_SAMPLE] + sample_idx;
                    int imap = sample_idx;
                    subset_rec = bcf_dup(var);
                    bcf_subset(hdr, subset_rec, 1, &imap);

                    m_collectors[sample->key]->processVariant(m_subset_hdrs[sample_names[sample_idx]], subset_rec);

                    bcf_destroy(subset_rec);
                }

            }

            virtual void appendJsonImpl(json_t * jsonRootObj) override {
                for(auto& sample : m_collectors) {
                    json_t *j_sample = json_object();
                    sample.second->appendJson(j_sample);
                    json_object_set_new(jsonRootObj, sample.first.c_str(), j_sample);
                }
            }

        public:
            BySampleStratifier(bcf_hdr_t* hdr) : AbstractStatCollector() { 
                auto sample_iter = hdr->id[BCF_DT_SAMPLE];
                auto sample_iter_end = hdr->id[BCF_DT_SAMPLE] + hdr->n[BCF_DT_SAMPLE];

                sample_count = hdr->n[BCF_DT_SAMPLE];
                sample_names = new char*[hdr->n[BCF_DT_SAMPLE]];

                for(int sample_idx = 0; sample_idx < hdr->n[BCF_DT_SAMPLE]; sample_idx++) {
                    int imap = sample_idx;
                    sample_names[sample_idx] = strdup(hdr->id[BCF_DT_SAMPLE][sample_idx].key);
                    std::cerr<<"Sample ["<<sample_names[sample_idx]<<"] seen; creating collector"<<std::endl;
                    m_collectors[sample_names[sample_idx]] = new CollectorT();
                    m_subset_hdrs[sample_names[sample_idx]] = bcf_hdr_subset(hdr, 1, hdr->samples + sample_idx, &imap);
                }

                
            }
            virtual ~BySampleStratifier() {
                // release collectors for each sample
                for(auto& sample : m_collectors) delete sample.second;

                // release bcf_hdr_t objects for each sample
                for(auto& sample : m_subset_hdrs) bcf_hdr_destroy(sample.second);

                // release sample names and pointer array
                for(int i=0; i<sample_count; i++) free(sample_names[i]);
                delete [] sample_names;
            }

        private:
            std::map<std::string, CollectorT*> m_collectors;
            int sample_count;
            char ** sample_names;
            std::map<std::string, bcf_hdr_t*> m_subset_hdrs;
    };
}

#endif
