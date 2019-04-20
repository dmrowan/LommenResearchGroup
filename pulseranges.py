#!/usr/bin/env python
import argparse
import numpy as np
import LCClass

desc="""
Find primary and interpulse phase ranges with user
defined offpeak phase ranges
"""

def find_phase_ranges(evt, off1, off2, nsigma):
    lc = LCClass.LightCurve(evt)
    lc.generate()
    cutofftup = lc.peak_cutoff(off1, off2, nsigma=nsigma)
    print(f"Min phase primary: {cutofftup.min_phase}")
    print(f"Max phase primary: {cutofftup.max_phase}")
    print(f"Min phase interpulse: {cutofftup.min_phase_ip}")
    print(f"Max phase interpulse: {cutofftup.max_phase_ip}")
    return cutofftup


if __name__ == '__main__':
    parser = arparse.ArgumentParser(description=desc)
    parser.add_argument("--evt", help="Event file", type=str, 
                        required=True)
    parser.add_argument("--lower", help="Lower off peak phase", 
                        type=float, required=True)
    parser.add_argument("--upper", help="Upper off peak phase",
                        type=float, required=True)
    parser.add_argument("--sigma", help="Sigma limits for cutoff", 
                        type=float, default=3.0)
    args= parser.parse_args()
    find_phase_range(args.evt, args.lower, args.upper, args.sigma)
