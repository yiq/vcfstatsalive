#include "AbstractStatCollector.h"
#include "BasicStatsCollector.h"
#include "ByGenotypeStratifier.h"
#include "BySampleStratifier.h"

using namespace VcfStatsAlive;

class StatsCollector : public BasicStatsCollector {
    public:
        StatsCollector() : BasicStatsCollector(1, 200, false) {}
};

void printStatsJansson(AbstractStatCollector* rootStatCollector) {

	// Create the root object that contains everything
	json_t * j_root = json_object();

	// Let the root object of the collector tree create Json
	rootStatCollector->appendJson(j_root);

	// Dump the json
    std::cout<<json_dumps(j_root, JSON_COMPACT | JSON_ENSURE_ASCII | JSON_PRESERVE_ORDER)<<";"<<std::endl;

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

    //ByGenotypeStratifier<StatsCollector> strat;

    bcf_hdr_t* hdr = bcf_hdr_read(fp);

    BySampleStratifier<ByGenotypeStratifier<StatsCollector>> strat(hdr);

	bcf1_t* line = bcf_init();
    
	while(bcf_read(fp, hdr, line) == 0) {
        strat.processVariant(hdr, line);
    }

    printStatsJansson(&strat);

    bcf_destroy1(line);
    bcf_hdr_destroy(hdr);

    hts_close(fp);
}
