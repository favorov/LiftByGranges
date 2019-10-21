#!/usr/bin/env python
desc="""Convert BAM to BigWig."""

import os, sys, pysam
import pybedtools
from datetime import datetime

def bam_to_bigwig_edit(bam, genome, output, scale=False):
    """
    Modified from bam2bigwig.py from https://github.com/lpryszcz/bin so that it accepts
    unsorted BAM files for conversion to BigWig.
    
    Given a BAM file `bam` and assembly `genome`, create a bigWig file scaled
    such that the values represent scaled reads -- that is, reads per million
    mapped reads. (Disable this scaling step with scale=False; in this case values will
    indicate number of reads)

    Assumes that `bedGraphToBigWig` from UCSC tools is installed; see
    http://genome.ucsc.edu/goldenPath/help/bigWig.html for more details on the
    format.
    """
    genome_file = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    kwargs = dict(bg=True, split=True, g=genome_file)
    if scale:
        readcount = mapped_read_count(bam)
        _scale = 1 / (readcount / 1e6)
        kwargs['scale'] = _scale
    x = pybedtools.BedTool(bam).genome_coverage(**kwargs)
    y = x.sort()
    cmds = [
        'bedGraphToBigWig',
        y.fn,
        genome_file,
        output]
    os.system(' '.join(cmds))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--bam", required=True,
                        help="BAM file")
    parser.add_argument("-g", "--genome", required=True,
                        help="genome FASTA file")
    parser.add_argument("-o", "--output", required=True,
                        help="output stream [stdout]")
    parser.add_argument("--scaling", default=False,  action="store_false",
                        help="disable RPM scaling")
    o = parser.parse_args()
    sys.stderr.write("Options: %s\n"%str(o))
        
    bam_to_bigwig_edit(o.bam, o.genome, o.output, o.scaling)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed! \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
