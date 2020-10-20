from pysuffixarray.core import SuffixArray
from collections import Counter
import os


# 1/a is fraction of rows kept in suffix array sample, 1/b is the fraction of rows checkpointed
class FMIndex():
    def __init__(self, genome, a, b):
        genome = genome.lower()
        genome = genome.replace("$", "")
        a = abs(int(a))
        b = abs(int(b))
        self.sa = SuffixArray(genome).suffix_array() # creates suffix array
        # use suffix array to form bwt and suffix array sample
        self.bwt = self.bwt(genome)
        self.a = a
        self.sa_sample = self.sa_sample()
        del self.sa # now that we are done using the full suffix array, delete it from memory
        self.empty_character_counts = self.create_character_counts() # get the number of times each character occurs in a dictionary as an alternative to the first row of the bwt matrix
        self.b = b
        self.checkpoints = self.checkpoints() # set up checkpoints for how many times each character has occured at different locations in the bwt, also creates first occurence dict
        self.first_occs = self.create_first_occs()
        
    def bwt(self, genome): # get the bwt of the string from the suffix array
        genome = genome.lower()
        bwt = []
        for suffix_index in self.sa:
            if suffix_index == 0: 
                bwt.append('$')
            else: 
                bwt.append(genome[suffix_index-1])
        return ''.join(bwt)
    
    def sa_sample(self): # take 1/a values from the suffix array to the sample suffix array
        sa_sample = []
        for i, suffix_index in enumerate(self.sa):
            if i % self.a == 0:
                sa_sample.append(suffix_index)
        return sa_sample
        
    def create_character_counts(self): # set all character counts to 0
        unique_characters = sorted(list(set(self.bwt)))
        character_counts = dict()
        for char in unique_characters:
            character_counts[char] = 0
        return character_counts
    
    def checkpoints(self): # create checkpoints of character counts throughout bwt
        checkpoints = []
        self.character_counts = dict(self.empty_character_counts)
        for i in range(len(self.bwt) + 1):
            if i != 0: # when i = 0, no characters are bring counted
                self.character_counts[self.bwt[i-1]] += 1 # this also updates the characters counts, which give us all of the information for the first row
            if i % self.b == 0: # every 1/b characters will have a checkpoint
                checkpoints.append(dict(self.character_counts))
        return checkpoints
    
    def create_first_occs(self):
        position = 0
        first_occs = dict(self.empty_character_counts)
        for key, value in self.character_counts.items():
            first_occs[key] = position
            position += value
        return first_occs
    
    def get_checkpoint(self, pos):
        remainder = pos % self.b
        count_dict = dict(self.empty_character_counts)
        if remainder == 0:
            return Counter(self.checkpoints[pos//self.b])
        elif remainder < self.b/2 or (len(self.bwt)-1)//self.b == pos//self.b:
            for i in range(remainder):
                count_dict[self.bwt[pos-i-1]] += 1
            final_checkpoint = Counter(self.checkpoints[pos//self.b])
            final_checkpoint.update(Counter(count_dict))
            return final_checkpoint
        else:
            for i in range(self.b - remainder):
                 count_dict[self.bwt[pos+i]] += 1
            final_checkpoint = Counter(self.checkpoints[pos//self.b + 1])
            final_checkpoint.subtract(Counter(count_dict))
            return final_checkpoint
    
    
    def read_fastq(self, fastq):
        if not os.path.exists(fastq):
            return "The input fastq file does not exist. Please try again."
        
        ks = ['name', 'sequence', 'quality']
        
        with open(fastq, 'r') as file:
            sequences = []
            lines = []
            count = 0
            for line in file:
                count += 1
                if count != 3:
                    lines.append(line.rstrip())
                if len(lines) == 3:
                    record = {k:v for k,v in zip(ks, lines)}
                    sequences.append(record)
                    lines = []
                    count = 0
        return sequences
    
    
    def EXACTMATCH(self, fastq):
        
        fastq = self.read_fastq(fastq)
        if fastq == "The input fastq file does not exist. Please try again.":
        	return fastq
        for sequence in fastq:
            read = sequence['sequence']
        if len(read) == 0:
            return "Your read must be length 1 or greater."
        if len(read) > len(self.bwt):
            return "The length of the read should not be greater than the length of the genome."
        read = read.lower()
        try:
            i = self.first_occs[read[-1]]
        except KeyError:
            return "All of the characters in your read must also be present in the genome."
        j = self.first_occs[read[-1]] + self.character_counts[read[-1]]
        #mismatches = [0 for i in range(len(self.bwt))]
        for char in read[::-1][1:]:
            i_checkpoint = self.get_checkpoint(i)
            j_checkpoint = self.get_checkpoint(j)
            first_char_rank = i_checkpoint[char]
            char_count = j_checkpoint
            char_count.subtract(i_checkpoint)
            char_count = char_count[char]
            try:
                i = self.first_occs[char] + first_char_rank
            except KeyError:
                return "All of the characters in your read must also be present in the genome."
            j = i + char_count
            if i == j:
                return "0 matches found."
        matched_positions = []
        for index, char in enumerate(self.bwt[i:j]):
            pos = i + index
            count = 0
            while pos % self.a != 0:
                count += 1
                rank_dict = self.get_checkpoint(pos)
                pos = self.first_occs[char] + rank_dict[char]
                char = self.bwt[pos]
            matched_positions.append((self.sa_sample[pos//self.a] + count) % len(self.bwt))
            
        if j - i == 1:
            return "1 match found at position {}!".format(matched_positions[0])
        else:
            return "{} matches found at positions {}!".format(j-i, matched_positions)            
        
    
         
    
fm = FMIndex("ACGCTCTTCAG",1 , 1)
print(fm.EXACTMATCH("test.fq"))