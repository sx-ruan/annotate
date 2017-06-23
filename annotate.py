# coding=utf-8

"""
A simple Python program for annotating VCF files

Required python packages: requests, numpy
Usage: python annotate.py out_path [vcf_path]

The output file is tab-delimited and has 9 columns:
chrom: Chromosome
pos: Reference position
ref: Reference sequence
alt: Alternate (variant) sequence
variant_type: Variant type composed of one or two items. The first item comes 
    from the VCF file and is either snp, mnp, ins, del or complex. The optional 
    second item is a Sequence Ontology (SO) term that describes the consequence 
    of the variant, if the variant is in the ExAC database. If several SO terms 
    are possible, only the most severe one is reported.
seq_depth: Depth of sequence coverage at the site of variation
num_variant_reads: Number of reads supporting the variant
pct_variant_reads: Percentage of reads supporting the variant versus those 
    supporting the reference
exac_allele_freq: Allele frequency of variant from Broad Institute ExAC Project

Note:
If the VCF file records multiple alternate alleles at the same position, all of 
the variants will be recorded (but separately) in the annotation file.
"""

import csv
import requests
import warnings
import sys
import numpy as np

# Broad Institute ExAC API URL
EXAC_URL = 'http://exac.hms.harvard.edu/rest/variant/'

# Order of severity of Sequence Ontology (SO) terms (more severe to less severe)
# as estimated by Ensembl
# Source: http://www.ensembl.org/info/genome/variation/predicted_data.html
ENSEMBL_SO_SEVERITY_ORDER = (
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant'
)


class ExAC(dict):
    """
    Class for retrieving and parsing ExAC variant information
    """
    def __init__(self, variant_id):
        """
        :param variant_id: a string that identifies the variant: CHROMOSOME-POSITION-REFERENCE-VARIANT
        """
        r = requests.get(EXAC_URL + variant_id)
        super(ExAC, self).__init__(r.json())

    @classmethod
    def severity(cls, so_term):
        """
        Get the severity ranking of an SO term
        :param so_term: SO term
        :return: severity ranking
        """
        try:
            return ENSEMBL_SO_SEVERITY_ORDER.index(so_term)
        except ValueError:
            warnings.warn('Unexpected Sequence Ontology term: %s' % so_term)
            # If the SO term is not in ENSEMBL_SO_SEVERITY_ORDER, its severity
            # ranking is assumed to be +inf (least severe)
            return np.inf

    def allele_freq(self):
        """
        If the variant is in the ExAC database, return its allele frequency, 
        otherwise return NaN
        :return: ExAC allele frequency
        """
        return float(self['variant'].get('allele_freq', np.nan))

    def consequence(self):
        """
        If the variant is in the ExAC database, return its most severe SO term, 
        otherwise return an empty string
        :return: SO term or an empty string
        """
        d = self['consequence']
        return min(d.keys(), key=self.severity) if d else ''


class VCFReader(object):
    """
    Class for reading and parsing VCF v4.1 files
    """
    def __init__(self, path):
        """
        :param path: path of a VCF file
        """
        self.file = open(path)
        self.reader = csv.DictReader(self.parse(self.file), dialect='excel-tab')

    @classmethod
    def parse(cls, f):
        """
        Skip the meta-information lines in a VCF file
        :param f: an iterable
        """
        for row in f:
            if row.startswith('##'):
                continue
            elif row.startswith('#'):
                yield row.lstrip('#')  # preprocess the header row
            else:
                yield row

    class Info(dict):
        """
        Class for parsing the INFO column in a VCF file
        """
        def __init__(self, string):
            """
            :param string: an entry in the INFO column
            """
            super(VCFReader.Info, self).__init__()

            for field in string.split(';'):
                data = field.split('=')

                if len(data) == 1:  # handle INFO fields of the 'Flag' type
                    key, = data
                    val = True
                elif len(data) == 2:  # handle INFO fields of other types
                    key, val = data
                else:
                    raise ValueError

                self[key] = val


def main(out_path, vcf_path='Challenge_data.vcf'):
    """
    Generate the variant annotation file
    :param out_path: path of the output file
    :param vcf_path: path of the input VCF v4.1 file
    """
    vcf = VCFReader(vcf_path)  # read and parse the VCF file

    with open(out_path, 'wb') as f:
        writer = csv.writer(f, dialect='excel-tab')
        header = ('chrom', 'pos', 'ref', 'alt', 'variant_type', 'seq_depth',
                  'num_variant_reads', 'pct_variant_reads', 'exac_allele_freq')
        writer.writerow(header)  # write the header of the annotation file

        for row in vcf.reader:
            chrom = row['CHROM']
            pos = row['POS']
            ref = row['REF']
            alts = row['ALT'].split(',')  # list of variants

            info = VCFReader.Info(row['INFO'])  # parse the INFO entry
            info_dp = int(info['DP'])  # total read depth at the locus
            info_ro = float(info['RO'])  # reference allele read count
            info_aos = map(int, info['AO'].split(','))  # variant read counts

            # types of variants, either snp, mnp, ins, del, or complex
            info_types = info['TYPE'].split(',')

            assert len(alts) == len(info_aos) == len(info_types)

            for i, alt in enumerate(alts):
                variant_id = '-'.join((chrom, pos, ref, alt))  # get variant id
                exac = ExAC(variant_id)  # retrieve and parse ExAC variant data
                allele_freq = exac.allele_freq()
                consequence = exac.consequence()  # get the most severe SO term

                variant_type = ('%s, %s' % (info_types[i], consequence)
                                if consequence else info_types[i])

                # percentage of reads supporting the variant versus those
                # supporting the reference allele
                pct = info_aos[i] / (info_ro + info_aos[i])

                txts = (chrom, pos, ref, alt, variant_type, info_dp,
                        info_aos[i], '%#.4G' % pct, '%#.4G' % allele_freq)
                writer.writerow(txts)

    vcf.file.close()


if __name__ == '__main__':
    warnings.resetwarnings()
    warnings.simplefilter('always', UserWarning)
    main(*sys.argv[1:])
