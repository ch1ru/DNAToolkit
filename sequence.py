import collections
from models.codon import codons
from models.nucleotides import Nucleotides, Compliments

class Sequence:

    def __init__(self, seq: str):
        self.seq = self.validateSeq(seq)
        if not self.seq:
            raise Exception("Invalid sequence")

    @staticmethod
    def validateSeq(seq):
        
        tmp_seq = seq.upper()
        for n in tmp_seq:
            if n not in Nucleotides:
                return False
        return tmp_seq

    @staticmethod
    def countNucFrequency(seq, n: str = None):
        if n:
            return dict(collections.Counter(seq)).get(n)
        else:
            return dict(collections.Counter(seq))

    @staticmethod
    def transcribe(seq: str):
        seq = seq.upper()
        rna = seq.replace('T', 'U')
        return rna
    
    @staticmethod
    def reverse_compliment(seq):
        
        dna = ""
        for n in seq:
            dna += Compliments[n]
        return dna[::-1]
    
    @staticmethod
    def gc_calculation(seq):
        return round((seq.count('C') + seq.count('G')) / len(seq) * 100)

    @staticmethod
    def gc_content_substr(seq, k=20):
        res = []
        for i in range(0, len(seq) - k + 1, k):
            substr = seq[i:i + k]
            res.append(Sequence.gc_calculation(substr))
        return res

    @staticmethod
    def translate_seq(seq, init_pos=0):
        return [codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]
    
    @staticmethod
    def codon_usage(seq, aminoacid):
        tmp = []
        for i in range(0, len(seq) - 2, 3):
            if codons[seq[i: i + 3]] == aminoacid:
                tmp.append(seq[i:i + 3])
        
        freqDict = dict(collections.Counter(tmp))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict
    
    @staticmethod
    def gen_reading_frames(seq):
        frames = []
        frames.append(Sequence.translate_seq(seq, 0))        
        frames.append(Sequence.translate_seq(seq, 1))
        frames.append(Sequence.translate_seq(seq, 2))
        frames.append(Sequence.translate_seq(Sequence.reverse_compliment(seq), 0))
        frames.append(Sequence.translate_seq(Sequence.reverse_compliment(seq), 1))
        frames.append(Sequence.translate_seq(Sequence.reverse_compliment(seq), 2))
        return frames
    
    @staticmethod
    def proteins_from_rf(aa_seq):
        curr_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                #STOP codon
                if curr_prot:
                    for p in curr_prot:
                        proteins.append(p)
                    curr_prot = []
            else:
                #START codon
                if aa == 'M':
                    curr_prot.append("")
                    for i in range(len(curr_prot)):
                        curr_prot[i] += aa
        return proteins

