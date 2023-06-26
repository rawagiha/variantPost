class Variant(object):
    def __init__(self, chrom, pos, ref, alt, reference, window=0):

        self.chrom = chrom
        self.pos = pos  # 1-based
        self.ref = ref
        self.alt = alt
        self.reference = reference

        self.window = 200 if window <= 200 else window

        # chrom name chrom end check
        self._chrom = self.chrom

        self.reference_len = reference.get_reference_length(self._chrom)

        # 1-based
        self.unspliced_local_reference_start = max(0, pos - self.window) + 1
        self.unspliced_local_reference_end = min(pos + self.window, self.reference_len)

        #self.unspliced_local_reference = reference.fetch(
        #    chrom,
        #    self.unspliced_local_reference_start - 1,
        #    self.unspliced_local_reference_end,
        #)

    @property
    def is_leftaligned(self):
        if self.ref[-1].upper() != self.alt[-1].upper():
            return True
        elif "N" in self.ref.upper() or "N" in self.alt.upper():
            return True
        
        return False

    @property
    def is_normalized(self):
        if self.is_leftaligned:
            if len(self.ref) > 1 and len(self.alt) and (self.ref[0].upper() == self.alt[0].upper()):
                return False
            else:
                return True
        else:
            return False

    
    def normalize(self, inplace=False):

        if inplace:
            i = self
        else:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference)

        condition_1 = i.ref[-1].upper() == i.alt[-1].upper() != "N"
        lhs = i.reference.fetch(
            i.chrom, max(0, i.pos - 1 - int(self.window * 1.5)), i.pos - 1
        )[::-1]
        
        n = 0
        lhs_len = len(lhs)
        while condition_1 and n < lhs_len:

            left_base = lhs[n]

            i.ref = left_base + i.ref[:-1]
            i.alt = left_base + i.alt[:-1]
            i.pos -= 1

            condition_1 = i.ref[-1].upper() == i.alt[-1].upper() != "N"

            n += 1

        condition_2 = i.ref[0].upper() == i.alt[0].upper()
        condition_3 = len(i.ref) > 1 and len(i.alt) > 1
        while condition_2 and condition_3:
            i.ref = i.ref[1:]
            i.alt = i.alt[1:]
            i.pos += 1
            condition_2 = i.ref[0].upper() == i.alt[0].upper()
            condition_3 = len(i.ref) > 1 and len(i.alt) > 1

        if inplace:
            return None
        else:
            return i
