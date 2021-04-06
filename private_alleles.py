import os 
import argparse
from collections import defaultdict

import vcf
import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(
            namespace, self.dest, os.path.abspath(os.path.expanduser(values))
        )


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Get private alleles""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--base",
        required=True,
        type=is_file,
        action=FullPaths,
        help="""The base vcf.gz file.""",
    )
    parser.add_argument(
        "--alternate",
        required=True,
        type=is_file,
        action=FullPaths,
        help="""The alternate vcf.gz file""",
    )
    parser.add_argument(
        "--base-name",
        type=str,
        default=None,
        help="""A different name for the base population""",
    )
    parser.add_argument(
        "--alternate-name",
        type=str,
        default=None,
        help="""A different name for the alternate population""",
    )

    return parser.parse_args()


def make_data_dictionary(vcf_file):
    data = defaultdict(lambda: defaultdict(dict))
    for record in vcf.Reader(filename=vcf_file, compressed=True):
        #{chrom:
        #    {pos:
        #            'alleles':{'REF':{},'ALT':{}}
        #            'genotypes':set()
        #    }
        calls = {'REF':record.REF, 'ALT':record.ALT}
        genotypes = set()
        alleles = set()
        who = defaultdict(set)
        for sample in record.samples:
            if '.' in sample['GT']:
                pass
            else:
                cleaned_gt = sample['GT'].replace("|","/")
                genotypes.add(cleaned_gt)
                gt_list = cleaned_gt.split("/")
                alleles.update(gt_list)
                for allele in gt_list:
                    who[allele].add(sample.sample)
        data[record.CHROM][record.POS] = {
            'calls':calls,
            'genotypes':genotypes,
            'alleles':alleles,
            'who':who
        }
    return data

def output_comparison(base, alt, basename=None, altname=None, verbose=True):
    count = 0
    print("chrom;pos;calls;pvt_alleles;base_geno({});alt_geno({});base_who_count;base_who".format(
        basename,
        altname
    ))
    for chrom, data in base.items():
        for pos, info in data.items():
            base_alleles = info['alleles']
            base_genotypes = info['genotypes']
            alt_alleles = alt[chrom][pos]['alleles']
            alt_genotypes = alt[chrom][pos]['genotypes']
            # genotypes in base not in alt
            diff = base_alleles.difference(alt_alleles)
            if diff != set() and verbose:
                sum_of_who = []
                #pdb.set_trace()
                for allele in list(diff):
                    sum_of_who.extend(info['who'][allele])
                print("{};{};{};{};{};{};{};{}".format(
                    chrom,
                    pos,
                    info['calls'],
                    diff,
                    base_genotypes,
                    alt_genotypes,
                    len(sum_of_who),
                    sum_of_who
                ))
                count += 1
    #print("Total count ({} v. {})= {}".format(basename, altname, count))


def main():
    args = get_args()
    base = make_data_dictionary(args.base)
    alt = make_data_dictionary(args.alternate)
    if not args.base_name:
        base_name = os.path.basename(args.base)
    else:
        base_name = args.basename
    if not args.alternate_name:
        alternate_name = os.path.basename(args.alternate)
    else:
        alternate_name = args.alternate_name
    count = 0
    output_comparison(base, alt, base_name, alternate_name)
    print("\n\n")
    output_comparison(alt, base, alternate_name, base_name)

if __name__ == '__main__':
    main()