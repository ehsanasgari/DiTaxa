'''
Created on Feb 5, 2016

@author: ehsan
'''
from Bio.Blast import NCBIWWW
from taxonomic_assignment.blastresult import BLASTResults

class BLASTER(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''

    @staticmethod
    def dnaSeqBLAST(sequence,E_VALUE_THRESH = 0.001):
        '''
            Return BLAST results on a DNA sequence
        '''
        return BLASTResults(NCBIWWW.qblast("blastn", "nt", sequence),E_VALUE_THRESH)

    @staticmethod
    def blast16S(sequence,E_VALUE_THRESH = 0.00001):
        '''
            Return BLAST results on a DNA sequence
        '''
        return BLASTResults(NCBIWWW.qblast("blastn", 'MG1655', sequence),E_VALUE_THRESH)

    @staticmethod
    def protSeqBLAST(sequence,E_VALUE_THRESH = 0.001):
        '''
            Return BLAST results on a protein sequence
        '''
        return BLASTResults(NCBIWWW.qblast("blastp", "nr", sequence, hitlist_size=1000), E_VALUE_THRESH)

    @staticmethod
    def dnaFastaBLAST(fastaAddress, hitlist_size=1000,E_VALUE_THRESH = 0.002, batch_size=20):
        '''
            Return the results on a batch of DNA sequences
        '''
        with open(fastaAddress) as myfile:
            lines=myfile.readlines()
            BR=BLASTResults()
            for i in range(0,len(lines),batch_size):    
                BR.addResult(NCBIWWW.qblast("blastn", "nt", ''.join(lines[i:(i+batch_size)]), hitlist_size=hitlist_size),E_VALUE_THRESH)
        return BR

    @staticmethod
    def protFastaBLAST(fastaAddress, hitlist_size=10,E_VALUE_THRESH = 0.001, batch_size=20):
        '''
            Return the results on a batch of protein sequences
        '''  
        with open(fastaAddress) as myfile:
            lines=myfile.readlines()
            BR=BLASTResults()
            for i in range(0,len(lines),batch_size):    
                BR.addResult(NCBIWWW.qblast("blastp", "nr", ''.join(lines[i:(i+batch_size)]), hitlist_size=hitlist_size),E_VALUE_THRESH)
        return BR

    @staticmethod
    def mergeMultipleProteinHits(listOfSequenceDic, header=[]):
        '''
            This function return list of list of sequences in different species
        '''
        species_lists = [list(sequenceDic.keys()) for sequenceDic in listOfSequenceDic]
        common_species=list(set(species_lists[0]).intersection(*species_lists))
        sequencelist=[]
        for sequenceDic in listOfSequenceDic:
            sequencelist.append([sequenceDic[key].replace("-","") for key in common_species])
        #species names
        common_species=[k.replace(" ","_") for k in common_species]
        return common_species,sequencelist
    
    @staticmethod
    def mergeAligments(algn1,algn2, header=[]):
        '''
            This function return list of list of sequences in different species
        '''
        common_species=list(set(algn1.keys()).intersection(set(algn2.keys())))
        finalAlignments=[]
        for species in common_species:
            finalAlignments.append(algn1[species]+algn2[species])
        return finalAlignments,common_species

