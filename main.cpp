#include <iostream>
#include <getopt.h>

#include <string>
#include <fstream>

#include <memory>

#include <execinfo.h>

#include "BasicStatsCollector.h"

using namespace std;
using namespace VcfStatsAlive;
using namespace htslib;

static struct option getopt_options[] =
{
	/* These options set a flag. */
	{"update-rate",		required_argument,	0, 'u'},
	{"first-update",	required_argument,	0, 'f'},
	{"qual-lower-val",	optional_argument,	0, 'q'},
	{"qual-upper-val",	optional_argument,	0, 'Q'},
	{"log-scale-af",	optional_argument,	0, 'l'},
	{0, 0, 0, 0}
};

static unsigned int updateRate;
static unsigned int firstUpdateRate;
static int qualHistLowerVal;
static int qualHistUpperVal;

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

void printStatsJansson(AbstractStatCollector* rootStatCollector);

int main(int argc, char* argv[]) {

	signal(SIGSEGV, handler);   // install our handler


	string filename;
	updateRate = 1000;
	firstUpdateRate = 0;
	qualHistLowerVal = 1;
	qualHistUpperVal = 200;
	bool logScaleAF = false;

	int option_index = 0;

	int ch;
	while((ch = getopt_long (argc, argv, "f:u:q:Q:l", getopt_options, &option_index)) != -1) {
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

	vcf::VariantCallFile vcfFile;
	htsFile *fpVcf;

	if (argc == 0) {
		vcfFile.open(std::cin);
		fpVcf = hts_open("_", "r");
	}
	else {
		filename = *argv;
		vcfFile.open(filename);
		printf("Opening file '%s'\n", filename.c_str());
		fpVcf = hts_open(filename.c_str(), "r");
	}

	if(!vcfFile.is_open()) {
		std::cerr<<"Unable to open vcf file / stream"<<std::endl;
		exit(1);
	}

	BasicStatsCollector *bsc = new BasicStatsCollector(qualHistLowerVal, qualHistUpperVal, logScaleAF);


	unsigned long totalVariants = 0;
	vcf::Variant var(vcfFile);

	while(vcfFile.is_open() && !vcfFile.done()) {

		vcfFile.getNextVariant(var);
		bsc->processVariant(var);
		totalVariants++;

		if((totalVariants > 0 && totalVariants % updateRate == 0) ||
				(firstUpdateRate > 0 && totalVariants >= firstUpdateRate)) {

			printStatsJansson(bsc);

			// disable first update after it has been fired.
			if(firstUpdateRate > 0) firstUpdateRate = 0;
		}
	}

	printStatsJansson(bsc);

	delete bsc;

	return 0;
}

void printStatsJansson(AbstractStatCollector* rootStatCollector) {

	// Create the root object that contains everything
	json_t * j_root = json_object();

	// Let the root object of the collector tree create Json
	rootStatCollector->appendJson(j_root);

	// Dump the json
	cout<<json_dumps(j_root, JSON_COMPACT | JSON_ENSURE_ASCII | JSON_PRESERVE_ORDER)<<";"<<endl;

	json_decref(j_root);
}
