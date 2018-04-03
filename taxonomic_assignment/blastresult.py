'''
Created on Feb 8, 2016

@author: ehsaneddin
'''
from Bio.Blast import NCBIXML

class BLASTResults(object):
    '''
    classdocs
    '''


    def __init__(self,blast_result=None,E_VALUE_THRESH = 0.001):
        '''
        Constructor
        '''
        # [length of alignment, score of alignment, query_start_pos, query_ene_pos, entropy of sequence_query ]
        self.result_rows=[]
        # mapping from GID to the description
        self.speciesGID2name=dict()
        # mapping from GID to row in the results
        self.speciesGID2row=dict()
        # mapping from row ID to GID
        self.rowidx2GID=[]
        # Add the results
        if blast_result!=None:
            self.addResult(blast_result, E_VALUE_THRESH)
        
    def addResult(self,blast_result,E_VALUE_THRESH = 0.001):
        '''
        This method adds new results
        '''
        # Parsing the results  
        blast_records_list=[];
        blast_records_list.append(NCBIXML.parse(blast_result))
       
        for blast_records in blast_records_list:
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            # if the result has the certain threshold of e_value add it
                            self.result_rows.append([(hsp.align_length - hsp.gaps), hsp.score , hsp.query_start, hsp.query_end, hsp.query, hsp.match, hsp.sbjct])
                            # species_strain assignment
                            if alignment.hit_id in self.speciesGID2row:
                                self.speciesGID2row[alignment.hit_id].append(len(self.result_rows)-1)
                            else:
                                self.speciesGID2row[alignment.hit_id]=[len(self.result_rows)-1]
                                self.speciesGID2name[alignment.hit_id]=alignment.hit_def
                            
                            self.rowidx2GID.append(alignment.hit_id)
    
    def getSpeciesSequences(self):
        '''
            This function extract the sequences and the species from the BLAST results
        '''
        matchedPairs=dict()
        for k,v in self.speciesGID2row.items():
            vname=self.speciesGID2name[k]
            species=vname[vname.find("[")+1:vname.find("]")]
            species=species.replace(" ","_")
            row=self.result_rows[v[0]]
            sequence=row[7]
            matchedPairs[species]=sequence
        return matchedPairs
