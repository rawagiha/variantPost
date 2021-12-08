
class Variant(object):
    def __init__(self, chrom, pos, ref, alt, reference, first_creation=True, window=0):

        self.chrom = chrom
        self.pos = pos    # 1-based
        self.ref = ref
        self.alt = alt
        self.reference = reference 
         
        self.window = 200 if window <= 200 else window

        #chrom name chrom end check 
        self._chrom = self.chrom
        
        self.reference_len = reference.get_reference_length(self._chrom) 
        
        # 0-based
        self.unspliced_local_reference_start = max(0, pos - self.window) 
        self.unspliced_local_reference_end = min(pos + self.window, self.reference_len)
        
        self.unspliced_local_reference = reference.fetch(
            chrom, self.unspliced_local_reference_start, self.unspliced_local_reference_end
        )

