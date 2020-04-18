import collections
import logging
import sys

import vcf

contig_spec = collections.namedtuple('Contig', 'id,length')
calldata_spec = collections.namedtuple('CallData', 'GT')
CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']

def convert_hg19_vcf_to_grch37_vcf(input_vcf_file, output_vcf_file):
    """
    Liftover lifts from GRCh38 -> hg19. hg19 contigs have 'chr' prefix,
    that needs to be removed for GRCH37. Also filtering to only acceptable chrs
    autosomes, X and Y. No alternate contigs are considered
    :param input_vcf_file: hg19 vcf file liftedover from GRCh38
    :param output_vcf_file: GRCh37 contigs vcf file
    """
    chroms_to_keep = ['chr' + chrom for chrom in CHROMS]
    grch37_contigs = collections.OrderedDict()
    with open(input_vcf_file) as in_fp, open(output_vcf_file, 'w') as out_fp:
        reader = vcf.Reader(in_fp)
        # update contigs to remove 'chr' prefix
        for contig in reader.contigs:
            if contig in chroms_to_keep:
                value = contig_spec(reader.contigs[contig].id.replace('chr', ''), reader.contigs[contig].length)
                grch37_contigs.update({contig.replace('chr', ''): value})

        reader.contigs = grch37_contigs
        writer = vcf.Writer(out_fp, reader, lineterminator='\n')
        for record in reader:
            if record.CHROM in chroms_to_keep:
                # update chrom in the records
                record.CHROM = record.CHROM.replace('chr', '')
                writer.write_record(record)


MISMATCH_SITES = [
    {
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 21012603,
            'end': 21012604,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 21235475,
            'end': 21235476,
            'base': 'T'
        }
    },
    {
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 7563750,
            'end': 7563751,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 7563983,
            'end': 7563984,
            'base': 'T'
        }
    },
    {
        '38_coordinates': {
            'chrom': 'chr15',
            'start': 48515440,
            'end': 48515441,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr15',
            'start': 48807637,
            'end': 48807638,
            'base': 'C'
        }
    },
    {
        '38_coordinates': {
            'chrom': 'chr19',
            'start': 55154216,
            'end': 55154217,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr19',
            'start': 55665584,
            'end': 55665585,
            'base': 'A'
        }
    },
]


def record_overlaps_mismatch_sites(record):
    for site in MISMATCH_SITES:
        if record.CHROM == site['38_coordinates']['chrom']:
            if record.POS <= site['38_coordinates']['start'] <= record.end:
                return site
    return False


def update_grch38_ref_to_grch37_for_record_if_needed(record, mismatched_site=None):
    """
    If record overlaps mismatched sites, update ref, alt and genotype accordingly
    Some assumptions:
    - records have only one sample (true for color data, one file per sample)
    - will not try to fix record, if genotype is malformed or missing
    - will not try to fix record, if more than 2 alleles
    """
    if not mismatched_site:
        mismatched_site = record_overlaps_mismatch_sites(record)
    if not mismatched_site:
        return record

    # get genotype
    # code only for one sample per vcf file. This is how color files are
    assert (len(record.samples) == 1)
    sample_call = record.samples[0]
    if 'GT' not in record.FORMAT:
        # missing genotype issue not handled
        raise ValueError(f'GT not in record {record}')
    gt_indices = sample_call.gt_alleles
    expected_gts = set(str(i) for i in range(len(record.ALT) + 1))
    if not set(gt_indices).issubset(expected_gts):
        # unknown genotype issue not handled
        raise ValueError(f'Unknown genotypes for record {record}, genotype {gt_indices}')

    ref = list(record.REF)
    pos = mismatched_site['38_coordinates']['start'] - record.POS

    if mismatched_site['38_coordinates']['base'] == ref[pos]:
        ref[pos] = mismatched_site['37_coordinates']['base']
    updated_ref = ''.join(ref)
    orig_allele_seqs = []
    # alleles from original record
    for allele in sample_call.gt_alleles:
        # only add reference, if it is in genotype, het variant
        if allele == '0':
            orig_allele_seqs.append(record.REF)
        else:
            orig_allele_seqs.append(record.ALT[int(allele) - 1])
    updated_alts = []
    updated_gt = []
    for allele_seq in orig_allele_seqs:
        # original allele matches updated ref, genotype needs '0' added
        if allele_seq == updated_ref:
            updated_gt.append(0)
        else:
            if allele_seq not in updated_alts:
                updated_alts.append(allele_seq)
            updated_gt.append(updated_alts.index(allele_seq) + 1)
    updated_gt = '/'.join(map(str, sorted(updated_gt)))

    if len(updated_alts) > 2:
        # case we will not try to handle, write original record as is
        raise ValueError(f'Updated record has more than 2 alts {updated_alts}, original record {record}')

    if len(updated_alts) == 0:
        # return None if the variant was not a variant in 37
        record = None
    else:
        # add genotype data to record. Loses all other fields
        record.FORMAT = 'GT'
        record.samples[0].data = calldata_spec(updated_gt)

        record.REF = updated_ref
        record.ALT = []
        for alt in updated_alts:
            record.ALT.append(vcf.model._Substitution(alt))
    return record


def convert_grch38_ref_mismatch_sites_to_grch37(input_vcf_file, output_vcf_basename):
    """
    For ACMG59 reportable range there are 4 sites that have
    reference mismatch between GRCh37 and GRCh38
    All ref and alts in variants overlapping these sites
    will need to be updated to 37 reference
    2 output vcf files will be created, one with variants overlapping mismatchsites
    and one with all other variants with original record
    """
    logger = logging.getLogger(__name__)
    ref_mimatch_vcf = f'{output_vcf_basename}_ref_mismatch.vcf'
    output_vcf_file = f'{output_vcf_basename}.vcf'
    reader = vcf.Reader(open(input_vcf_file))
    records = list(reader)
    mismatched_site_overlap = {}
    with open(output_vcf_file, 'w') as out_fp, open(ref_mimatch_vcf, 'w') as updated_fp:
        writer = vcf.Writer(out_fp, reader, lineterminator='\n')
        updated_writer = vcf.Writer(updated_fp, reader, lineterminator='\n')
        for record in records:
            mismatched_site = record_overlaps_mismatch_sites(record)
            if mismatched_site:
                mismatched_site_overlap[str(mismatched_site)] = True
                try:
                    new_record = update_grch38_ref_to_grch37_for_record_if_needed(record, mismatched_site)
                except ValueError as e:
                    logger.info(e)
                    new_record = record
                if new_record:
                    updated_writer.write_record(new_record)
            else:
                writer.write_record(record)
        # if there are no overlapping variants in mismatched sites,
        # create a homozygous variant matching 37 as ref and 38 as alt
        record.FORMAT = 'GT'
        record.INFO = []
        record.samples[0].data = calldata_spec('1/1')
        for site in MISMATCH_SITES:
            if str(site) not in mismatched_site_overlap.keys():
                record.CHROM = site['38_coordinates']['chrom']
                record.POS = site['38_coordinates']['start']
                record.REF = site['37_coordinates']['base']
                record.ALT = [vcf.model._Substitution(site['38_coordinates']['base'])]
                updated_writer.write_record(record)


if __name__ == '__main__':
    function = getattr(sys.modules[__name__], sys.argv[1])
    function(sys.argv[2], sys.argv[3])
