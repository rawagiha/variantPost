import numpy as np
from difflib import SequenceMatcher
from collections import OrderedDict, Counter

from .variant import Variant

def phase(contig_dict, skips, target_pos):
    if len(contig_dict) == 0:
        return "aho"
    
    snvs, indels, actual_events = profile_non_refs(contig_dict, target_pos)

    if actual_events.is_null():
        raise ValueError("invalid contig")

    contig_dict, snvs, indels = crop_contig(contig_dict, skips, snvs, indels, actual_events.pos)
    
    if not snvs and not indels:
        print("simple, nothin to phase")
        return "simple, nothing to phase" 
    else:
        print("may be complex")
        return "may be complex"
    #for k, v in contig_dict.items():
    #    print(k, v[0], v[1]) 


class NonReferenceEvent(object):
    def __init__(self, pos, ref, alt, qual):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = qual
    
    def __eq__(self, other):
        pos_eq = (self.pos == other.pos)
        ref_eq = (self.ref == other.ref)
        alt_eq = (self.alt == other.alt)

        return all([pos_eq, ref_eq, alt_eq])

    def is_null(self):
        if self.ref == self.alt == "N":
            return True
        else:
            return False


def parse_contig_for_exons(contig_dict, skips):
    contig_start = list(contig_dict.items())[0][0]
    contig_end = list(contig_dict.items())[-1][0]
    
    exon_start = contig_start
    exon_end = contig_end
    exons = []
    for skip in skips:
         exon_end = skip[0] - 1
         exons.append((exon_start, exon_end))
         exon_start = skip[1] + 1
         exon_end = contig_end
    
    exons.append((exon_start, exon_end))
    
    return contig_start, contig_end, exons



def profile_non_refs(contig_dict, target_pos):
    actual_pos = -1
    actual_event = NonReferenceEvent(-1, "N", "N", 0)
    pos_diff = np.inf

    snvs, indels = [], []
    for pos, v in contig_dict.items():
        ref, alt, qual = v[0], v[1], v[2]

        if ref != alt:
            
            non_ref = NonReferenceEvent(pos, ref, alt, qual)
            
            if abs(pos - target_pos) < pos_diff:
                actual_pos = pos
                pos_diff = abs(pos - target_pos)
                actual_event = non_ref

            if len(ref) == len(alt):
                snvs.append(non_ref)
            else:
                indels.append(non_ref)

    if actual_event in snvs:
        snvs.remove(actual_event)

    if actual_event in indels:
        indels.remove(actual_event)
    
    return snvs, indels, actual_event    


def crop_contig(contig_dict, skips, snvs, indels, target_pos):
    
    lt_lim, rt_lim, exons = parse_contig_for_exons(contig_dict, skips)
    for exon in exons:
        if exon[0] <= target_pos <= exon[1]:
            lt_lim = exon[0] - 1
            rt_lim = exon[1] + 1
   
    lt_disqualified, rt_disqualified = [], []
    for pos, v in contig_dict.items():
        if "N" in v[0] or "N" in v[1]:
            if pos < target_pos:
                lt_disqualified.append(pos)
            else:
                rt_disqualified.append(pos)

    if lt_disqualified:
        lt_lim = max(lt_lim, max(lt_disqualified)) 

    if rt_disqualified:
        rt_lim = min(rt_lim, min(rt_disqualified))

    contig_dict = {pos : v for pos, v in contig_dict.items() if lt_lim < pos < rt_lim}
    snvs = [snv for snv in snvs if lt_lim < snv.pos < rt_lim]
    indels = [indel for indel in indels if lt_lim < indel.pos < rt_lim]

    return contig_dict, snvs, indels


