from .variant import Variant, to_minimal_repeat_unit, repeat_counter

revcomp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


class IndelTaxon(object):
    def __init__(self, ref, alt, lt_flank, rt_flank):
        self.class_89 = "NA"
        self.class_83 = "NA"

        ## pre-left align + pre-complex filter + pre-snv.mnv filter

        if len(ref) < len(alt):
            self._type = "Ins"
            self.seq = alt[1:]
        else:
            self._type = "Del"
            self.seq = ref[1:]

        self.len = len(self.seq)
        self.unit = ""
        self.ulen = -1
        if self._type == "Ins":
            # Do not count repeats in the inserted seq
            self.rep = 0
            # Count by seq
            self.rep += repeat_counter(self.seq, rt_flank)
        else:
            self.unit = to_minimal_repeat_unit(self.seq)
            self.ulen = len(self.unit)
            # Count repeats in the deleted seq
            self.rep = repeat_counter(self.unit, self.seq)
            # Count by unit
            self.rep += repeat_counter(self.unit, rt_flank)

        self.mlen = -1
        if self._type == "Del" and self.len > 1 and self.rep == 1:
            i = 0
            is_same = self.seq[i] == rt_flank[i]
            while is_same:
                i += 1
                is_same = self.seq[i] == rt_flank[i]

            if i > 0:
                self.mlen = i

        self.flank_5p = ""
        self.flank_3p = ""
        if self.len == 1:
            self.flank_5p = lt_flank[-1]

            i = 0
            is_same = self.seq == rt_flank[i]
            while is_same:
                i += 1
                is_same = self.seq == rt_flank[i]

            if i < len(rt_flank):
                self.flank_3p = rt_flank[i]

            else:
                self.flank_3p = "N"

            if self.seq in ("G", "A"):
                self.seq = revcomp[self.seq]
                tmp = self.flank_5p
                self.flank_5p = revcomp[self.flank_3p]

                self.flank_3p = revcomp[tmp]

        if self._type == "Ins":
            if self.len == 1:
                if self.seq == "C":
                    self.single_C_insertion_classify()
                else:
                    self.single_T_insertion_classify()
            else:
                self.longer_insertion_classify()
        else:
            if self.len == 1:
                if self.seq == "C":
                    self.single_C_deletion_classify()
                else:
                    self.single_T_deletion_classify()
            else:
                if self.mlen >= 1:
                    self.micro_homology_deletion_classify()
                else:
                    self.longer_non_microhomo_deletion_classify()

    def single_C_insertion_classify(self, tuck_in=True):
        rep_83 = min(5, self.rep)
        self.class_83 = "1:Ins:C:{}".format(rep_83)

        if 0 <= self.rep <= 3:
            self.class_89 = "Ins(C):R(0,3)"

            if self.rep == 0 and self.flank_5p == "A":
                if self.flank_3p == "A":
                    self.class_89 = "A[Ins(C):R0]A"
                elif self.flank_3p == "T":
                    self.class_89 = "A[Ins(C):R0]T"
        elif 4 <= self.rep <= 6:
            self.class_89 = "Ins(C):R(4,6)"
        elif 7 <= self.rep <= 9:
            self.class_89 = "Ins(C):R(7,9)"
        else:
            # treat +9
            if tuck_in:
                self.class_89 = "Ins(C):R(7,9)"
            else:
                pass

    def single_T_insertion_classify(self, tuck_in=True):
        rep_83 = min(5, self.rep)
        self.class_83 = "1:Ins:T:{}".format(rep_83)

        if 0 <= self.rep <= 4:
            rep_89 = "R(0,4)"
        elif 5 <= self.rep <= 7:
            rep_89 = "R(5,7)"
        elif 8 <= self.rep <= 9:
            rep_89 = "R(8,9)"
        else:
            if tuck_in:
                rep_89 = "R(8,9)"
            else:
                return

        self.class_89 = "{}[Ins(T):{}]{}".format(self.flank_5p, rep_89, self.flank_3p)

    def single_C_deletion_classify(self, tuck_in=True):
        rep_83 = min(5, self.rep - 1)
        self.class_83 = "1:Del:C:{}".format(rep_83)

        if 1 <= self.rep <= 3:
            if self.flank_3p in ("A", "T"):
                self.class_89 = "[Del(C):R{}]{}".format(self.rep, self.flank_3p)
            else:
                self.class_89 = "[Del(C):R(1,5)]G"
        elif 4 <= self.rep <= 5:
            if self.flank_3p in ("A", "T"):
                self.class_89 = "[Del(C):R(4,5)]{}".format(self.flank_3p)
            else:
                self.class_89 = "[Del(C):R(1,5)]G"
        elif 6 <= self.rep <= 9:
            self.class_89 = "Del(C):R(6,9)"
        else:
            if tuck_in:
                self.class_89 = "Del(C):R(6,9)"

    def single_T_deletion_classify(self, tuck_in=True):
        rep_83 = min(5, self.rep - 1)
        self.class_83 = "1:Del:T:{}".format(rep_83)

        if 1 <= self.rep <= 4:
            self.class_89 = "{}[Del(T):R(1,4)]{}".format(self.flank_5p, self.flank_3p)
        elif 5 <= self.rep <= 7:
            self.class_89 = "{}[Del(T):R(5,7)]{}".format(self.flank_5p, self.flank_3p)
        elif 8 <= self.rep <= 9:
            self.class_89 = "{}[Del(T):R(8,9)]{}".format(self.flank_5p, self.flank_3p)
        else:
            if tuck_in:
                self.class_89 = "{}[Del(T):R(8,9)]{}".format(
                    self.flank_5p, self.flank_3p
                )

    def longer_insertion_classify(self, tuck_in=True):
        rep_83 = min(5, self.rep)
        len_83 = min(5, self.len)

        self.class_83 = "{}:Ins:R:{}".format(len_83, rep_83)

        if self.rep == 0:
            if 2 <= self.len <= 4:
                self.class_89 = "Ins(2,4):R0"
            else:
                self.class_89 = "Ins(5,):R0"
        elif self.rep == 1:
            if 2 <= self.len <= 4:
                self.class_89 = "Ins(2,4):R1"
            else:
                self.class_89 = "Ins(5,):R1"
        elif 2 <= self.rep <= 4:
            self.class_89 = "Ins(2,):R(2,4)"
        elif 5 <= self.rep <= 9:
            self.class_89 = "Ins(2,):R(5,9)"
        else:
            if tuck_in:
                self.class_89 = "Ins(2,):R(5,9)"

    def longer_non_microhomo_deletion_classify(self, tuck_in=True):
        rep_83 = min(5, self.rep - 1)
        len_83 = min(5, self.len)

        self.class_83 = "{}:Del:R:{}".format(len_83, rep_83)

        if self.rep == 1:
            if 2 <= self.len <= 4:
                self.class_89 = "Del(2,4):R1"
            else:
                self.class_89 = "Del(5,):R1"
        elif 2 <= self.rep <= 4:
            if 1 <= self.ulen <= 2:
                if self.len in (2, 3, 4, 6, 8):
                    self.class_89 = "Del(2,8):U(1,2):R(2,4)"
            else:
                if self.rep == 2:
                    self.class_89 = "Del(3,):U(3,):R2"
                else:
                    self.class_89 = "Del(3,):U(3,):R(3,9)"
        elif 5 <= self.rep <= 9:
            if 1 <= self.ulen <= 2:
                self.class_89 = "Del(2,):U(1,2):R(5,9)"
            else:
                self.class_89 = "Del(3,):U(3,):R(3,9)"
        else:
            if tuck_in:
                if 1 <= self.ulen <= 2:
                    self.class_89 = "Del(2,):U(1,2):R(5,9)"
                else:
                    self.class_89 = "Del(3,):U(3,):R(3,9)"

    def micro_homology_deletion_classify(self):
        len_83 = min(5, self.len)
        mh_83 = min(4, self.mlen)

        self.class_83 = "{}:Del:M:{}".format(len_83, mh_83)

        if self.mlen == 1:
            if 2 <= self.len <= 5:
                self.class_89 = "Del(2,5):M1"
            else:
                self.class_89 = "Del(6,):M1"
        elif self.mlen == 2:
            if 3 <= self.len <= 5:
                self.class_89 = "Del(3,5):M2"
            else:
                self.class_89 = "Del(6,):M2"
        elif self.mlen == 3:
            if 4 <= self.len <= 5:
                self.class_89 = "Del(4,5):M(3,4)"
            else:
                self.class_89 = "Del(6,):M3"
        else:
            if 4 <= self.len <= 5:
                self.class_89 = "Del(4,5):M(3,4)"
            else:
                self.class_89 = "Del(6,):M(4,)"
