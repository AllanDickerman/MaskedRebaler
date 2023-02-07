#!/usr/bin/env python3
"""
Mask final assembly (set low-coverage regions to Ns) for testing mask_low_coverage_regions
"""
import sys

from rebaler.mask_low_coverage_regions import mask_low_coverage_regions
 
if __name__ == '__main__':

    assembly, reads = sys.argv[1:3]
    min_depth = 5
    if len(sys.argv) > 4:
        min_depth = int(sys.argv[4])
    threads = 4
    if len(sys.argv) > 4:
        threads = int(sys.argv[4])
    sys.stderr.write("assembly={}, reads={}, min_depth={}, threads={}\n".format(assembly, reads, min_depth, threads))

    mask_low_coverage_regions(assembly, reads, min_depth=min_depth, threads=threads)
