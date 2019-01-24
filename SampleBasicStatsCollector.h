#ifndef SAMPLEBASICSTATSCOLLECTOR_H
#define SAMPLEBASICSTATSCOLLECTOR_H

#pragma once

#include "BasicStatsCollector.h"

namespace VcfStatsAlive {
    class SampleBasicStatsCollector : public BasicStatsCollector {
        public:
            SampleBasicStatsCollector(int qualLower, int qualUpper, bool logScaleAF):
                BasicStatsCollector(qualLower, qualUpper, logScaleAF) {}

            virtual void processVariantImpl(bcf_hdr_t* hdr, bcf1_t* var) override {
                // increment total variant counter
                ++_stats[kTotalRecords];

                int32_t ngt, *gt_arr=nullptr, ngt_arr=0;
                ngt = bcf_get_genotypes(hdr, var, &gt_arr, &ngt_arr);
                if (ngt <= 0) return; // no genotype info present

                bool isSnp = bcf_is_snp(var);
                int refLength = strlen(var->d.allele[0]);

                for(int altIndex = 0; altIndex < ngt_arr; altIndex++) {
                    int gt = (gt_arr[altIndex] >> 1) - 1;
                    if(gt <= 0) continue; // either missing or reference

                    updateTsTvRatio(var, gt, isSnp);
                    updateMutationSpectrum(var, gt, isSnp);
                    updateVariantTypeDist(var, gt, refLength);
                }

                updateAlleleFreqHist(hdr, var);
                updateQualityDist(var->qual);

                free(gt_arr);
            }
    };
};
#endif
