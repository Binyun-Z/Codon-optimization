import numpy as np
from numpy import random
import scipy.stats
import pandas as pd
import collections
from random import randint
from Bio.SeqUtils import GC
from Bio import SeqIO
import re

# ------------------------ Variable definition (codon table) --------------------------- #

# Assign full codon table on which to operate codon selection based on user's input
# 分配完整的密码子表，根据用户输入操作密码子选择
codon_table_full = { 'A': ['GCT','GCC','GCA','GCG'],
                'C': ['TGT','TGC'],
                'D': ['GAT','GAC'],
                'E': ['GAA','GAG'],
                'F': ['TTT','TTC'],
                'G': ['GGT','GGC','GGA','GGG'],
                'H': ['CAT','CAC'],
                'I': ['ATT','ATC','ATA'],
                'K': ['AAA','AAG'],
                'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
                'M': ['ATG'],
                'N': ['AAT','AAC'],
                'P': ['CCT','CCC','CCA','CCG'],
                'Q': ['CAA','CAG'],
                'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
                'S': ['AGT','AGC','TCT','TCC','TCA','TCG'],
                'T': ['ACT','ACC','ACA','ACG'],
                'V': ['GTT','GTC','GTA','GTG'],
                'W': ['TGG'],
                'Y': ['TAT','TAC'],
                '*': ['TAG','TAA','TGA'] }  


class DG:
    def __init__(self, cfg = None,codon_table_full = codon_table_full):
        self.my_prot = cfg.seqs #m目的蛋白
        self.seqVars =cfg.seqVars # 生成的不同基因数量
        self.host = cfg.host   # 表达宿主
        self.rar_thr=cfg.rar_thr #输入 RSCU 值，低于该值的密码子将被丢弃
        self.gc_thr = cfg.gc_thr #输入 RSCU 值，低于该值不以 G/C 结尾的密码子将被丢弃：
        self.codon_table_full = codon_table_full
        self.host_path = self.host+'_release.csv'
        self.codon_table = self.generate_codon_table()
    def generate_codon_table(self):
            rscu = pd.read_csv(self.host_path)
    
            # Extract amino acid symbols (one-letter code) from .csv file
            symb = rscu['AmOneLet']
            # Initialise a dictionary with unique amino acid symbols as keys
            codon_table = {}
            for s in symb:
                if s not in codon_table:
                    codon_table[s] = []
                    
            # Extract parameters from .csv file
            cod = rscu['Codon']
            val = rscu['RSCU']
            gc3 = rscu['GC3']
            
            # Fill codon table with codons above thresholds
            for index in range(len(symb)):
                # Only consider codons if above 'rare threshold'
                if val[index] > self.rar_thr:
                    # If above threshold, take all cods with rscu above gc_thr or ones ending in g/c (regardless of val) 
                    if val[index] > self.gc_thr or gc3[index] == 'Y':
                        codon_table[symb[index]].append(cod[index])
                        
                                
            # If RSCU thresholds set too high, most highly expressed codons are added to codon_table
            
            # Create a tuple with (symb, cod, val), used below to take most highly expressed codon
            bund = [ (symb[i], cod[i], val[i]) for i in range(len(cod)) ]
            
            # Cycle through keys (amino acids) in codon_table
            for key in codon_table:
                # If amino acids without corresponding codons are found..
                if codon_table[key] == []:
                    if key != 'W' and key != 'M':
                        print('Using only most highly expressed codon for amino acid %s' % key)
                    # Add most highly expressed codon from codon usage table
                    pool = [ bund[x] for x in range(len(bund)) if bund[x][0] == key ]
                    c = max(pool, key=lambda item:item[2])[1]
                    codon_table[key].append(c)
                    
            return codon_table
    def RSCU(self,path):
        """Given a path (argv) to codon usage table of expression host, functions generates a dictionary
        with all RSCU values (return)."""
        
        # Read codon_usage table from data folder
        codon_usage = pd.read_csv(path)
        
        # Initialise dictionary where RSCU values will be stored {triplet:RSCU val}
        RSCU_dict = {}
        
        # Cycle through codon_usage, extract info and load onto dictionary
        for i in range(len(codon_usage['Codon'])):
            RSCU_dict[codon_usage['Codon'][i]] = codon_usage['RSCU'][i]
            
        return  RSCU_dict
    def relative_adaptiveness(self,codon_table_full, path):
        """Calculates the relative adaptiveness of each codon given a path to a codon usage table."""
        
        # Define RSCU_dict using RSCU function
        RSCU_dict = self.RSCU(path)
        
        # For each amino acid, compute the max(RSCU) of synonymous codons
        # Initialise a dictionary, keys are cods, vals max RSCU
        maxRSCU_dict = {}
        for aa in codon_table_full:
            RSCUvals = [ RSCU_dict[cod] for cod in codon_table_full[aa] ]
            for cod in codon_table_full[aa]:
                maxRSCU_dict[cod] = max(RSCUvals)
        
        # Now define a relative adaptiveness dict, keys are cods, vals are ratios RSCU_cod/RSCU_max
        relative_adaptiveness_dict = {}
        for cod in RSCU_dict:
            relative_adaptiveness_dict[cod] = RSCU_dict[cod]/maxRSCU_dict[cod]
            
        return relative_adaptiveness_dict
    def CAI_calculator(self,codon_table_full, path, seq):

        # Generate relative_adaptiveness dictionary
        relative_adaptiveness_dict = self.relative_adaptiveness(codon_table_full, path)

        if len(seq) % 3 != 0:
            raise ValueError('length of sequence is not a multiple of three')
        
        # Initialise a list where codon adaptiveness values will be stored     
        seq_adaptiveness = []

        # Traverse coding sequence and extract codons
        for i in range(0,len(seq),3):
            cod = seq[i:i+3]
            cod_adaptiveness = relative_adaptiveness_dict[cod]
            seq_adaptiveness.append(cod_adaptiveness)
        
        # Calculate the geometric mean for the list seq_adaptiveness (= CAI)
        CAI = scipy.stats.mstats.gmean(seq_adaptiveness)
        
        return CAI
    def find_identities(self,l):
        """
        Takes in a list and returns a dictionary with seqs as keys and positions of identical elements in list as values.
        argvs: l =  list, e.g. mat[:,x]
        """

        # the number of items in the list will be the number of unique types
        uniq = [item for item, count in collections.Counter(l).items()]
        
        # Initialise a dictionary that will hold the results
        identDict = {}
        for item in uniq:
            identDict[item] = [ x for x in range(len(l)) if l[x] == item ]
            
        return identDict 
    
    def alternate(self):
        """A function that produces a generator
        Usage: alternator = alternate() """ 

        while True:
            yield 0
            yield 1
    def lib_generator(self):
        """
        Function to generate set of diversified coding sequences. 
        生成一组多样化编码序列的函数。
        argvs: my_prot = amino acid sequence of the protein in question (string),
        seqVars = number of coding sequence to generate (integer),
        codon_table: codon usage table for the specified expression host (dict)
        """
        my_prot, seqVars, codon_table = self.my_prot,self.seqVars,self.codon_table
        ### --- Preliminary operations, e.g. variables' initiation
        # N is the number of codons, same as the number of aa's in our protein sequence
        # N 是密码子的数量，与我们蛋白质序列中 aa 的数量相同
        N = len(my_prot)
        
        # Initialise data structure to store CDS variants
        # 初始化数据结构以存储 CDS 变体
        # use np array with dtype = string of 3 chars, as each codon is stored separately
        # 使用 dtype = string of 3 chars 的 np array，因为每个密码子都是单独存储的
        mat = np.zeros((seqVars,N), dtype='S3') #s3表示长度为3的字符串
        
        # create a list (codVars) with all codon variants at each position
        codVars = []
        # Fill codVars with a loop
        for aa in my_prot:
            codVars.append(codon_table[aa])
            
        # Check that dimensions of codVars and recipient data structure mat are compatible
        if len(codVars) != mat.shape[1]:
            raise ValueError('list of codon vars and matrix incompatible')   
            
        ### --- 现在填充mat的前两个位置
        
        # Pos 1, where we ensure that ATG is the starting codon
        # Pos 1，我们确保 ATG 是起始密码子
        if codVars[0] != ['ATG']:
            raise ValueError('first codon is not ATG')
        mat[:,0] = codVars[0]

        # Pos 2, where we ensure available codon variants are spread across growing sequences
        # Pos 2，我们确保可用的密码子变体分布在不断增长的序列中
        cur = 0
        while cur < seqVars:
            # Trick to cycle through codon variants and assign them to mat
            # 循环密码子变体并将它们分配给mat
            codIndex = cur - (cur//len(codVars[1]))*len(codVars[1])
            mat[cur,1] = codVars[1][codIndex]
            cur += 1  
            
        ### --- Method for comparing sequences and sort them into different sets
        
        # cycle through codVars (= codon variants, list of lists) and fill mat 
        for k in range(2,len(codVars) ):
            
            # Initialise final and initial positions for joining codons
            # iniPos will change within the while loop
            finPos = k 
            iniPos = k - 1
            
            # Initialise data frame to hold identity dictionaries, could be a list of dicts
            idents = []
            
            
            #### From here we have an iterative logic block, whereby
            #### 1) identities of stretches with increasing length are evaluated;
            #### 2) for each position to fill (k), a list of dictionaries is returned
            
            # Create a switch to exit while loop when appropriate
            switch = 1
            
            while switch == 1:
            
                # Initialise joint Codons list
                jointCods = []
                
                # Create the set that will be submitted to find_identities() function
                for x in range(seqVars):
                    jointCods.append("".join([mat[x,p].decode() for p in range(iniPos,finPos)]))   
                    
                # Deploy find_identities() on jointCods
                identDict = self.find_identities(jointCods)         
                
                ## Add hypothetical case where there are no identities at all
                
                # Clause to break out of while loop: if there are no identities or if first
                # codon has been reached 
                
                ### Note: we do not want the last entry in idents to be singlets
                ### UNLESS, iniPos = k - 1, i.e. there no terminal identities
                
                if iniPos == k - 1 and max([len(x) for x in identDict.values()]) == 1:
                    ### We are allowed to have singlets and return
                    idents.append(identDict)
                    ### Exit while loop by turning switch off
                    switch = 0
                    
                ### We do not want to append the last comparison with singlets, so we do
                ### not append with the following condition, just turn off the switch    
                if max([len(x) for x in identDict.values()]) == 1 or iniPos == 0:
                    
                    ### Exit while loop by turning switch off
                    switch = 0
                    
                elif max([len(x) for x in identDict.values()]) > 1:
                    ### Add ident dict to idents
                    idents.append(identDict)
                    ### Do not turn switch off but move iniPos back instead
                    iniPos = iniPos - 1
                    
                else:
                    raise ValueError('problem with identDict evaluation or logic')
            
            
            ########### ----------------- Block 1 ends here ----------------------############
            
            # Case where only one codon is available: just fill all seqs with that
            if len(codVars[k]) == 1:
                cod = codVars[k][0] 
                for pos in range(seqVars):
                    mat[pos,k] = cod
                    
            else:
                ########### ----------------- Block 2 ----------------- ######################
                
                # Start with special condition, note list around list
                if len(idents) == 1:
                    outerDict_temp = {}
                    for key in idents[0]:
                        outerDict_temp[key] = [ idents[0][key] ]
                    
                if len(idents) > 1:
                    # Extract outer idents layer in dictionary
                    outerDict = idents[0]
                    outerDict_temp = {}
                    # Create a second dictionary with same keys as above, but empty value-lists
                    for item in outerDict:
                        outerDict_temp[item] = [] 
                        
                    # Sequence clustering
                    
                    # Cycle through all layers in idents (including 3-bp homology layer, i.e. outermost layer)
                    for cur in range(1, len(idents) + 1):
                        # print 'cur is:', cur
                        index = -cur
                        # Take sequences for which there are > 1 ID's
                        
                        longestHom = [ x for x in idents[index] if len(idents[index][x]) > 1 ]
                        # Cycle through items in longestHom   
                        for item in longestHom:
                    
                            # Define last triplet as key for outerDict_temp (storage Dict)
                            currentKey = item[-3:]
                        
                            # Define elements corresponding to item in longestHom
                            l = [ x for x in idents[index][item] ]
                        
                            # if outerDict_temp[currentKey] is empty, just append
                            if outerDict_temp[currentKey] == []:
                                outerDict_temp[currentKey].append(l)
                                
                            else:
                                
                        
                                # Two cases - 1) no item in l is present in outerDict_temp[currentKey]
                            
                                # all items in outerDict_temp[currentKey] as a set
                                all_items = set( x for sublist in outerDict_temp[currentKey] for x in sublist )
                            
                                # if l and all_items do not share elements
                                if set(l).isdisjoint(all_items):
                            
                                    # Append list l to outerDict_temp
                                    outerDict_temp[currentKey].append(l)
                            
                                # 2) they are not disjoint
                                else:
                                    
                                    # Cycle through sublists in dict
                                    for sublist in outerDict_temp[currentKey]:
                                        
                                        # if a sublist shares items with l
                                        if not set(l).isdisjoint(set(sublist)):
                                    
                                            # Find items NOT in common
                                            dif = set(l).difference(set(sublist))
                                        
                                            # Cycle through them
                                            for e in dif:
                                        
                                                # If e is not already in all_items
                                                if e not in all_items:
                                                
                                                    # Append element to sublist and update all_items
                                                    sublist.append(e)
                                                    all_items.add(e)
                                                    
                ## Special condition - if singlets are present throughout the depth of idents,
                ## these are not considered. Fix by collecting them at the end
                for key in outerDict_temp:
                    if outerDict_temp[key] == []:
                        # print 'FIXING BUG'
                        outerDict_temp[key] = [idents[0][key]]
                                        
                ####################### -------------------------- ###########################                
                    
                ########### Module 3 -  filling module ###################
                ### Now we work on outerDict_temp
                
                # Initialise generator/alternator
                #print 'alternator is being initialised'
                alternator = self.alternate() 
                #print 'alternator is on: ', alternator
                
            
                codList = [ x for x in codVars[k] ]    
                
                # Shuffle codList
                random.shuffle(codList)

                # Split the shuffled list into two
                x1 = codList[0:len(codList)//2]
                x2 = codList[len(codList)//2:len(codList)]
                # Then we put them together
                pool = [x1, x2]
            
                # Initialise state, only done on first iteration
                if k == 2:
                    state = alternator.__next__()
                    
                # Before starting work on this 'k', ensure state is on 0
                if state == 1:
                    state = alternator.__next__()
                
                # Cycle through items in dictionary, each set defined as workingList
                for item in outerDict_temp:
                    workingList = outerDict_temp[item]     ###### SET OF LISTS
                
                    # Cycle through sublists in workingList
                    for sub in workingList:
                        # On making this transition, we re-define x1 and x2
                    
                        for id in sub:                     ###### SINGLE ID      
                    
                            ### Within the same sub list, we can flick between x1 and x2 when they are used up
                            
                            if pool[state] == []:
                                ### Switch to codon sublist not in use
                                state = alternator.__next__()
                            
                                ### Reload previous state
                                pool[state - 1] = [ x for x in codVars[k] if x not in pool[state] ]
                                
                            ## Standard lines for choosing codon and appending to matrix
                            if pool[state] == []:
                                state = alternator.__next__()
                            cod = random.choice(pool[state])
                            mat[id,k] = cod
                        
                            pool[state].remove(cod)
                            
                        ## Here we re-define x1 and x2, based on state
                        ######### Note: now adopting a new strategy, previous cod in unitTest if nec.
                        if state == 0:
                            # x2 becomes = leftovers of x1 + x2, and state switches
                            x2 = pool[state] + pool[state - 1]
                            x1 = [ x for x in codVars[k] if x not in x2]
                            pool = [x1,x2]
                            # Now we switch the state
                            state = alternator.__next__()
                            # 
                        if state == 1:
                            # x2 is leftovers of x2, state remains the same
                            x2 = pool[state]
                            x1 = [ x for x in codVars[k] if x not in x2]
                            pool = [x1,x2]
                            #### note: remaining in state 1
            
                            
            # Alternator must be initialised "manually" in the special case where there is only one
            # codon option at k = 2.              
            if k == 2 and len(codVars[k]) == 1:
                alternator = self.alternate()
                state = alternator.__next__()
                
                    

        # Convert matrix mat to list of sequences
        mySeqs = []
        for ind in range(seqVars):
            s = ''
            for cod in mat[ind,:]:
                s = s + cod.decode()
            mySeqs.append(s)   
                
        return mySeqs
    def hamming_distance(self,s1,s2):
        """ Takes two strings as input and
        returns the number of mismatches as type integer"""
        if len(s1) != len(s2):
            raise ValueError("Undefined for sequences of unequal length")
        return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
    
    def hamming_matrix(self,pool):
        """Takes in a list of sequences (pool) and returns 
        a matrix with pairwise hamming distances"""
        
        # Make sure pool variable is a list
        if type(pool) != list:
            raise ValueError('Pool is not a list')
        
        # Array dimension
        dimLen = len(pool)
        
        # Generate array of zeros
        hammingMatrix = np.zeros((dimLen,dimLen))
        
        # Fill hammingMatrix with hamming distances - symmetrical matrix
        for i,ele_1 in enumerate(pool):
        
            for j,ele_2 in enumerate(pool):
            
                # Only fill in upper triangle
                if j >= i:
                    break
                # Calculate distance
                misMatches = self.hamming_distance(ele_1,ele_2)
                # Fill in array
                hammingMatrix[i,j] = misMatches
                # Same for symmetrical element
                hammingMatrix[j,i] = misMatches
                
        return hammingMatrix
    def hamming_stats(self,pool):
        """Generates a hammingMatrix with hamming_matrix()
        and calculates mean distance, minimum and maximum 
        divergence (in %)"""
        
        # Generate array with hamming distances
        divArray = self.hamming_matrix(pool)
        # Store len of divArray
        dimA = divArray.shape[0]
        
        # Extract upper triangular elements
        iu = np.triu_indices(dimA)
        upTri = [ ele for ele in divArray[iu] if ele != 0 ] # excludes diagonal values
        
        # Calculate mean of upTri and express it in percentage
        meanHam = np.mean(upTri)
        
        # Express mean in percentage - meanHam/len(sequence)
        meanHamPer = float( meanHam/(len(pool[0]))) * 100.0
        
        # Extract the smallest hamming distance
        minHam = np.min(upTri)
        
        # Express min in percentage
        minHamPer = float( minHam/(len(pool[0]))) * 100.0
        
        maxHam = np.max(upTri)
        
        # Express max in percentage - meanHam/len(sequence)
        maxHamPer = float( maxHam/(len(pool[0]))) * 100.0
        
        
        return (meanHamPer,minHamPer,maxHamPer)
    def longest_cont(self,s1,s2):
        """Works out the longest stretch of identical bases between two
        degenerate coding sequences"""
        if len(s1) != len(s2):
            raise ValueError("Undefined for sequences of unequal length")
        res = ''
        for ch1, ch2 in zip(s1,s2):
            if ch1 == ch2:
                res += '1'
            else:
                res += '0'        
        #print res   ### Need to subdivide into chunks of contiguous 1's (and 0's)
                
        return np.max([len(x) for x in re.compile("(1+1)*").findall(res)])
    def longest_cont_matrix(self,pool):

        """Takes in a list of sequences (pool) and returns 
        a matrix with longest stretch of identity"""
        
        # Make sure pool variable is a list
        if type(pool) != list:
            raise ValueError('Pool is not a list')
        
        # Array dimension
        dimLen = len(pool)
        
        # Generate array of zeros
        longestMatrix = np.zeros((dimLen,dimLen))
        
        # Fill longestMatrix with identity values - symmetrical matrix
        for i,ele_1 in enumerate(pool):
        
            for j,ele_2 in enumerate(pool):
            
                # Only fill in upper triangle
                if j >= i:
                    break
                # Calculate distance
                idents = self.longest_cont(ele_1,ele_2)
                # Fill in array
                longestMatrix[i,j] = idents
                # Same for symmetrical element
                longestMatrix[j,i] = idents
                
        return longestMatrix
    def abs_longest(self,pool):
        """takes in a list of sequences (pool) and returns the longest stretch of homology
        that can be found in the set"""
        
        # Call longest_cont_matrix() and take the upper triangular region
        mat = np.triu(self.longest_cont_matrix(pool))
        store_seq = []
        # dump non-zero elements into store_seq
        for row in mat:
            for ele in row:
                if ele != 0:
                    store_seq.append(ele)
        # return max
        return np.max(store_seq) 
    # --------------------------- Main functions' definitions: END ------------------------- #
    def run(self):
        mySeqs = self.lib_generator()
        # 
        # print mySeqs
        # 
        # Process data 
        # 1) Generate ID list
        IDlist = [ 'seq' + str(i) for i in range(len(mySeqs)) ]
        # 
        # Generate CAI values
        CAIlist = []
        for seq in mySeqs:
            CAIlist.append(self.CAI_calculator(codon_table_full = self.codon_table_full, path = self.host_path, seq = seq))
        #     
        ## Generate GC values
        GClist = []
        for seq in mySeqs:
            GClist.append(GC(seq))
        #     
        ## Create dataframe with output sequences and statistics
        output_data = {'seq ID': IDlist, 'Sequence': mySeqs, 'CAI': CAIlist, 'GC %': GClist}
        df = pd.DataFrame(output_data)
        print( df )
# 
        ## Communicate hamming distance stats and longest stretch of homology
        print( 'Stats: Mean, minimum and maximum hamming distances in the sequence set are (per cent):', self.hamming_stats(mySeqs))
        print ('Stats: Longest stretch of homology between any two sequences (in bp):', self.abs_longest(mySeqs))
        # 
        df.to_csv('codingSequenceVariants.csv', sep = ',')

class Parameter:
    def __init__(self,my_prot = '' ,seqVars = 10,host = 'CLIB',description='用户自定义参数',
                 rar_thr = 0.5,gc_thr = 0.8):
        self.my_prot = my_prot
        self.seqVars =seqVars # 生成的不同基因数量
        self.host = host   # 表达宿主
        # self.host_path = self.host+'_codon_usage.csv' #同义密码子使用说明
        self.rar_thr=rar_thr #输入 RSCU 值，低于该值的密码子将被丢弃
        self.gc_thr = gc_thr #输入 RSCU 值，低于该值不以 G/C 结尾的密码子将被丢弃：
        self.description = description
        self.seqs = self.get_seq()
    def get_seq(self):
        seqlist = []
        for seqercord in SeqIO.parse(self.my_prot,"fasta"):
            seqlist.append(str(seqercord.seq))
        return(seqlist[0])

    def __str__(self):
        return f'{self.name}={self.value}'
    
if __name__ == "__main__":
    import argparse
    # 创建解析器
    parser = argparse.ArgumentParser(description='用户自定义参数')
    # 添加参数
    parser.add_argument('-p', dest ='my_prot', default='test.fasta',help='氨基酸序列的fasta文件,例如aaa.fasta 或aaa.fa,此文件必须和GD在同一个文件夹下')
    parser.add_argument('-n' ,dest = 'seqVars', type=int,default=10,help='设计的不同基因序列数量')

    # 解析参数
    args = parser.parse_args()
    cfg = Parameter(args.my_prot,args.seqVars)
    instance = DG(cfg,codon_table_full)
    instance.run()