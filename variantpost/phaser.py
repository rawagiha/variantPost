import numpy as np
from difflib import SequenceMatcher
from collections import OrderedDict, Counter

from .variant import Variant

def phase(contig_dict, target_pos):
    ref_seq, contig_seq, snvs, indels, actual_event = preprocess(contig_dict, target_pos)
    
    print(ref_seq, contig_seq, "hey") 



class NonReferenceEvent(object):
    def __init__(self, pos, ref, alt, qual):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = qual

def preprocess(contig_dict, target_pos):
    actual_pos = -1
    actual_event = NonReferenceEvent(-1, "N", "N", 0)
    pos_diff = np.inf

    ref_seq, contig_seq = "", ""
    snvs, indels = [], []
    for pos, v in contig_dict.items():
        ref, alt, qual = v[0].decode("utf-8"), v[1].decode("utf-8"), v[2]
        
        ref_seq += ref
        contig_seq += alt

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

    return ref_seq, contig_seq, snvs, indels, actual_event    
