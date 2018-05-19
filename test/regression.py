import json
import os
import re
import shutil
import subprocess
import unittest


class IntegrationTests(unittest.TestCase):

    # Json keys enumerated in test methods are annotated as follows.
    # The type of the value accessed by the json key is indicated in the comment following the key.
    # If value is also a json key, the type of the key->value map is given.
    # Because all json keys are strings, a comment in () is given to state the type the string really
    # represents. Elipsis is given when all remaining types are the same.
    # Keys that access a json key with complex structure just have a value type of json.


    platinum_results_file = '/Users/timoc/source/frameshift/vcfstatsalive/platnum-results.json'
    bad_results_file = '/Users/timoc/source/frameshift/vcfstatsalive/bad-results.json'

    assets = {
            'platinum-exome.vcf.gz' : 'platinum-exome.json',
            }
            

    @classmethod
    def setUpClass(cls):
        print("Processing variants to generate test data")

        shutil.rmtree('output/')
        os.mkdir('output/')

        for k, v in IntegrationTests.assets.items():
            IntegrationTests._process_variants(k, v)

    #-----------------------------------------------------------------------------
    # Test methods
    #-----------------------------------------------------------------------------

    def test_top_level(self):
        IntegrationTests._run_test_on(
                self._test_top_level,
                IntegrationTests.platinum_results_file,
                IntegrationTests.bad_results_file)
  
    def test_allele_frequency(self):
        IntegrationTests._run_test_on(
                self._regress_allele_frequency,
                IntegrationTests.platinum_results_file,
                IntegrationTests.bad_results_file,
                key='af_hist')
    
    def test_mutation_spectrum(self):
        IntegrationTests._run_test_on(
                self._regress_mutation_spectrum,
                IntegrationTests.platinum_results_file,
                IntegrationTests.bad_results_file,
                key='mut_spec')

    def test_varian_types(self):
        IntegrationTests._run_test_on(
                self._regress_variant_types,
                IntegrationTests.platinum_results_file,
                IntegrationTests.bad_results_file,
                key='var_type')
    
    def test_quality_distribution(self):
        IntegrationTests._run_test_on(
                self._regress_quality_distribution,
                IntegrationTests.platinum_results_file,
                IntegrationTests.bad_results_file,
                key='qual_dist')

    def test_indel_size(self):
        IntegrationTests._run_test_on(
                self._regress_indel_size,
                IntegrationTests.platinum_results_file,
                IntegrationTests.bad_results_file,
                key='indel_size')

    #-----------------------------------------------------------------------------
    # Internal test implementation
    #-----------------------------------------------------------------------------

    def _test_top_level(self, expected_json, observed_json):

        top_level_keys = [
            'TotalRecords', # float
            'TsTvRatio',    # float
            'af_hist',      # json
            'mut_spec',     # json
            'var_type',     # json
            'qual_dist',    # json
            'indel_size'    # str(int) -> int
            ]

        self._validate_keys(top_level_keys, expected_json, observed_json)

        self.assertEqual(expected_json['TotalRecords'], observed_json['TotalRecords'])
        self.assertEqual(expected_json['TsTvRatio'],    observed_json['TsTvRatio'])

    def _regress_allele_frequency(self, expected_json, observed_json):
       
        af_hist_keys = [
            'usingLogScaleAF', # bool, str(int) -> int
            'afHistBins'       # str(int) -> int
            ]

        self._validate_keys(af_hist_keys, expected_json, observed_json)

        using_log_scale = observed_json['usingLogScaleAF']
        self.assertEqual(bool, type(using_log_scale))
        self.assertEqual(expected_json['usingLogScaleAF'], using_log_scale)

        observed_histogram = observed_json['afHistBins']
        expected_histogram = expected_json['afHistBins']
        self.assertEqual(len(expected_histogram), len(observed_histogram))
        for k, v in expected_histogram.items():
            self.assertTrue(k in observed_histogram)
            self.assertEqual(v, observed_histogram[k])

    def _regress_mutation_spectrum(self, expected_json, observed_json):
        
        mut_spec_keys = [
            'A', # str -> list(int)
            'C', # ...
            'G',
            'T'
            ]

        self._validate_keys(mut_spec_keys, expected_json, observed_json)

        for bp in mut_spec_keys:
            self.assertEqual(4, len(observed_json[bp]))
            for i in range(0, 4):
                self.assertEqual(int, type(observed_json[bp][i]))
                self.assertEqual(expected_json[bp][i], observed_json[bp][i])

    def _regress_variant_types(self, expected_json, observed_json):
    
        var_type_keys = [
            'SNP', # str -> int
            'INS', # ...
            'DEL',
            'OTHER'
            ]
    
        self._validate_keys(var_type_keys, expected_json, observed_json)

        for k in var_type_keys:
            self.assertEqual(int, type(observed_json[k]))
            self.assertEqual(expected_json[k], observed_json[k])

    def _regress_quality_distribution(self, expected_json, observed_json):
       
        qual_dist_keys = [
            'qualHistLowerBound', # int
            'qualHistUpperBound', # int
            'lowerBin',           # int
            'upperBin',           # int
            'regularBins',        # str(int) -> int
            ]

        self._validate_keys(qual_dist_keys, expected_json, observed_json)

        for int_key in ['qualHistLowerBound', 'qualHistUpperBound', 'lowerBin', 'upperBin']:
            self.assertEqual(int, type(observed_json[int_key]))
            self.assertEqual(expected_json[int_key], observed_json[int_key])

        bin_start = expected_json['qualHistLowerBound']
        bin_end = expected_json['qualHistUpperBound']

        expected_bins = expected_json['regularBins']
        observed_bins = observed_json['regularBins']

        qual_keys = ["{0}".format(i) for i in range(bin_start - 1, bin_end - bin_start + 1)]

        self._validate_keys(qual_keys, expected_bins, observed_bins)

        for q in qual_keys:
            self.assertEqual(int, type(observed_bins[q]))
            self.assertEqual(expected_bins[q], observed_bins[q])

    def _regress_indel_size(self, expected_json, observed_json):

        # The keys are data-dependent in this case, so simply use the 
        # keys from the expected json
        sizes = expected_json.keys()
        self._validate_keys(sizes, expected_json, observed_json)

        for s in sizes:
            self.assertEqual(int, type(observed_json[s]))
            self.assertEqual(expected_json[s], observed_json[s])

    #-----------------------------------------------------------------------------
    # Helpers
    #-----------------------------------------------------------------------------

    @staticmethod
    def _run_test_on(method, expected_json_file, observed_json_file, key=None):
        for k, v in IntegrationTests.assets.items():

            expected_json = IntegrationTests._get_json('data/' + v)
            observed_json = IntegrationTests._get_json('output/' + v)

            if key is not None:
                method(expected_json[key], observed_json[key])
            else:
                method(expected_json, observed_json)

    @staticmethod
    def _get_json(json_file):
        json_data = None
        with open(json_file) as fp:
            json_data = json.load(fp)
    
        if json_data is None:
            raise ValueError("Issue loading json file")

        return json_data

    @staticmethod
    def _process_variants(vcf_file, output_file):
        proc = subprocess.Popen(['../vcfstatsalive', 'data/' + vcf_file], stdout=subprocess.PIPE)
        last_line = None

        for l in proc.stdout:
            last_line = l.decode()

        out, err = proc.communicate()
        
        with open('output/' + output_file, 'w') as fh:
            last_line = re.sub(';$', '', last_line)
            fh.write(last_line)
        fh.close()

    def _validate_keys(self, keys, expected_json, observed_json):
        self.assertEqual(len(expected_json), len(observed_json))
        self.assertEqual(len(keys), len(observed_json))

        for k in keys:
            if k not in expected_json:
                raise ValueError("Expected json missing known key '{0}'".format(k))
            self.assertTrue(k in observed_json, "Observed json missing key '{0}'".format(k))

