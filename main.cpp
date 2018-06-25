#include <iostream>
#include <getopt.h>

#include <string>
#include <fstream>

#include <memory>

#include "BasicStatsCollector.h"

using namespace std;
using namespace VcfStatsAlive;

static struct option getopt_options[] =
{
	/* These options set a flag. */
	{"update-rate",		required_argument,	0, 'u'},
	{"first-update",	required_argument,	0, 'f'},
	{"qual-lower-val",	optional_argument,	0, 'q'},
	{"qual-upper-val",	optional_argument,	0, 'Q'},
	{"log-scale-af",	optional_argument,	0, 'l'},
	{"subset-samples",	optional_argument,	0, 'S'},
	{"batch",			optional_argument,	0, 'b'},
	{0, 0, 0, 0}
};

static unsigned int updateRate;
static unsigned int firstUpdateRate;
static int qualHistLowerVal;
static int qualHistUpperVal;
static bool subsetSamples = false;
static bool batch = false;

void printStatsJansson(AbstractStatCollector* rootStatCollector, ostream& outstream);

bool bcf_gt_is_ref(int gt1, int gt2) {
	int gt = bcf_alleles2gt(bcf_gt_allele(gt1), bcf_gt_allele(gt2));
	return ( gt == 0 );
}

int *bcf_is_refs(bcf_hdr_t* hdr, bcf1_t* var)
{

	int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr);
	int *arr = new int[nsmpl];
	int32_t *gt_arr = NULL, ngt_arr = 0;

	ngt = bcf_get_genotypes(hdr, var, &gt_arr, &ngt_arr);

	int max_ploidy = ngt/nsmpl;
	for (i=0; i<nsmpl; i++)
	{
		int32_t *ptr = gt_arr + i*max_ploidy;
		for (j=0; j<max_ploidy; j+=2)
		{
			// if true, the sample has smaller ploidy
			if ( ptr[j]==bcf_int32_vector_end ) break;

			// missing allele
			if ( bcf_gt_is_missing(ptr[j]) ) continue;

			arr[i] = bcf_gt_is_ref(ptr[0], ptr[1]);
		}
	}

	free(gt_arr);

   	return arr;
}

int main(int argc, char* argv[]) {

	string filename;
	updateRate = 1000;
	firstUpdateRate = 0;
	qualHistLowerVal = 1;
	qualHistUpperVal = 200;
	bool logScaleAF = false;
	bool batch = true;

	int option_index = 0;

	int ch;
	while((ch = getopt_long (argc, argv, "f:u:q:Q:l:b:S", getopt_options, &option_index)) != -1) {
		switch(ch) {
			case 0:
				break;
			case 'u':
				updateRate = strtol(optarg, NULL, 10);
				break;
			case 'f':
				firstUpdateRate = strtol(optarg, NULL, 10);
				break;
			case 'q':
				qualHistLowerVal = strtol(optarg, NULL, 10);
				if(qualHistLowerVal < 0) {
					cerr<<"Invalid quality histogram lowerbound value "<<qualHistLowerVal<<endl;
					exit(1);
				}
				break;
			case 'Q':
				qualHistUpperVal = strtol(optarg, NULL, 10);
				if(qualHistUpperVal < 0) {
					cerr<<"Invalid quality histogram upperbound value "<<qualHistUpperVal<<endl;
					exit(1);
				}
				break;
			case 'l':
				logScaleAF = true;
				break;
			case 'S':
				subsetSamples = true;
				break;
			case 'b':
				batch = true;
				break;
			default:
				break;
		}
	}

	if(qualHistUpperVal < qualHistLowerVal) {
		cerr<<"Quality histogram upperbound is lower than lowerbound"<<endl;
		exit(1);
	}

	argc -= optind;
	argv += optind;

	htsFile *fp;

	if (argc == 0) {
		fp = hts_open("-", "r");
	}
	else {
		filename = *argv;
		fp = hts_open(filename.c_str(), "r");
	}

	if (fp == NULL) {
		std::cerr<<"Unable to open vcf file / stream"<<std::endl;
		exit(1);
	}

	unsigned long totalVariants = 0;

	bcf_hdr_t* hdr = bcf_hdr_read(fp);
	bcf1_t* line = bcf_init();
	int nsmpl = bcf_hdr_nsamples(hdr);

	BasicStatsCollector **bsc_arr = (BasicStatsCollector **)malloc(sizeof(BasicStatsCollector) * nsmpl);
	BasicStatsCollector *bsc;
	if (subsetSamples) {

	// 	// bsc = new BasicStatsCollector[nsmpl];
		for (int k=0; k < nsmpl; k++) {
			bsc_arr[k] = new BasicStatsCollector(qualHistLowerVal, qualHistUpperVal, logScaleAF);
		}
	} else {
		bsc = new BasicStatsCollector(qualHistLowerVal, qualHistUpperVal, logScaleAF);
	}

	while(bcf_read(fp, hdr, line) == 0) {

		// Unpack alternates and info block
		if (bcf_unpack(line, BCF_UN_ALL) != 0) {
			std::cerr<<"Error unpacking"<<std::endl;
		}


		if (subsetSamples) {
			int *isRefs = bcf_is_refs(hdr, line);
			for (int k=0; k < nsmpl; k++) {
				if (!isRefs[k]) {
					bsc_arr[k]->processVariant(hdr, line);
				}
			}

		} else {
			bsc->processVariant(hdr, line);
		}


		totalVariants++;

		if( (!subsetSamples || !batch) && ((totalVariants > 0 && totalVariants % updateRate == 0) ||
				(firstUpdateRate > 0 && totalVariants >= firstUpdateRate))) {

			printStatsJansson(bsc, cout);

			// disable first update after it has been fired.
			if(firstUpdateRate > 0) firstUpdateRate = 0;
		}
	}

	if(subsetSamples) {
		for (int k=0; k<nsmpl; k++) {
			ofstream myfile;
			char *sampleName = strcat(hdr->samples[k], ".json");
  			myfile.open (sampleName);
			printStatsJansson(bsc_arr[k], myfile);
			myfile.close();
		}
	} else {
		printStatsJansson(bsc, cout);
	}


	delete bsc;
	free(bsc_arr);

	return 0;
}

void printStatsJansson(AbstractStatCollector* rootStatCollector, ostream& outstream) {

	// Create the root object that contains everything
	json_t * j_root = json_object();

	// Let the root object of the collector tree create Json
	rootStatCollector->appendJson(j_root);

	// Dump the json
	outstream<<json_dumps(j_root, JSON_COMPACT | JSON_ENSURE_ASCII | JSON_PRESERVE_ORDER)<<";"<<endl;

	json_decref(j_root);
}
