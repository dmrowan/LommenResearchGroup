#!/usr/bin/env python
import argparse
import numpy as np
import LCClass

def find_phase_ranges(evt, off1, off2, nsigma):
    lc = LCClass.LightCurve(evt)
    lc.generate()
    cutofftup = lc.peak_cutoff(off1, off2, nsigma=nsigma)
    print(f"Min phase primary: {cutofftup.min_phase}")
    print(f"Max phase primary: {cutofftup.max_phase}")
    print(f"Min phase interpulse: {cutofftup.min_phase_ip}")
    print(f"Max phase interpulse: {cutofftup.max_phase_ip}")
    return cutofftup


