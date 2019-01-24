#include "AbstractStatCollector.h"
#include "SampleBasicStatsCollector.h"
#include "ByGenotypeStratifier.h"
#include "BySampleStratifier.h"

#include <csignal>

using namespace VcfStatsAlive;

class StatsCollector : public SampleBasicStatsCollector {
    public:
        StatsCollector() : SampleBasicStatsCollector(1, 200, false) {}
};


size_t line_count = 0;

void progressSignalHandler(int signum) {
    std::cerr<<line_count <<" vcf lines processed."<<std::endl;
}

void printStatsJansson(AbstractStatCollector* rootStatCollector) {

	// Create the root object that contains everything
	json_t * j_root = json_object();

	// Let the root object of the collector tree create Json
	rootStatCollector->appendJson(j_root);

	// Dump the json
    std::cout<<json_dumps(j_root, JSON_COMPACT | JSON_ENSURE_ASCII | JSON_PRESERVE_ORDER)<<std::endl;

	json_decref(j_root);
}

int main(int argc, const char** argv) {
    argc--; argv++;

    htsFile *fp;

    if (argc == 0) fp = hts_open("-", "r");
    else fp = hts_open(*argv, "r");


    if (fp == NULL) {
        if (argc == 0) std::cerr<<"Unable to open vcf stream"<<std::endl;
        else std::cerr<<"Unable to open vcf file "<<*argv<<std::endl;
        exit(1);
    };

    bcf_hdr_t* hdr = bcf_hdr_read(fp);

    BySampleStratifier<ByGenotypeStratifier<StatsCollector>> strat(hdr);

	bcf1_t* line = bcf_init();

    signal(SIGUSR1, progressSignalHandler);
    
	while(bcf_read(fp, hdr, line) == 0) {
        strat.processVariant(hdr, line);
        line_count++;
    }

    printStatsJansson(&strat);

    bcf_destroy1(line);
    bcf_hdr_destroy(hdr);

    hts_close(fp);
}
