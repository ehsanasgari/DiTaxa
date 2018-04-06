import sys
sys.path.append('../')
from utility.file_utility import FileUtility
from make_representations.cpe_apply import CPE
from make_representations.cpe_efficient import train_cpe
from multiprocessing import Pool
import tqdm
from sklearn.feature_extraction.text import TfidfVectorizer


def segmentseq(id_seq):
    global CPE_applier
    idx,seq=id_seq
    global CPE_applier
    return idx, CPE_applier.segement(seq)


train=FileUtility.load_list('~/datasets/plasmit/plasmid_training.txt')
test=FileUtility.load_list('~/datasets/plasmit/plasmid_testing.txt')

negative_class=FileUtility.read_fasta_sequences('/net/refdata/ncbi_2005/PROCESSED_NCBI/genomes/bacteria/bacteria_chromosomes_only.fasta')

f=open('../temp_db/plasmid_cpe','r')
global CPE_applier
CPE_applier=CPE(f,separator='', merge_size=305000)

all_sequences={idx:seq for idx,seq in enumerate(train+test+negative_class)}




pool = Pool(processes=20)
segmentation=dict()
for idx, val in tqdm.tqdm(pool.imap_unordered(segmentseq, list(all_sequences.items()), chunksize=5),
                       total=len(all_sequences)):
    segmentation[idx]=val

idxs=list(segmentation.keys())
idxs.sort()
corpus=[segmentation[x] for x in idxs]
tfm = TfidfVectorizer(use_idf=False, analyzer='word', tokenizer=str.split, ngram_range=(1,1), norm=None, stop_words=[], lowercase=False, binary=False)
tf_vec = tfm.fit_transform(corpus)
feature_names = tfm.get_feature_names()

FileUtility.save_sparse_csr('all_genomes',tf_vec)
FileUtility.save_list('all_genomes_features',feature_names)
FileUtility.save_list('all_genomes_lens',[str(x) for x in [len(train),len(test),len(negative_class)]])

tfm = TfidfVectorizer(use_idf=True, analyzer='word', tokenizer=str.split, ngram_range=(1,1), norm=None, stop_words=[], lowercase=False, binary=False)
tf_vec = tfm.fit_transform(corpus)
feature_names = tfm.get_feature_names()

FileUtility.save_sparse_csr('all_genomes_idf',tf_vec)
FileUtility.save_list('all_genomes_features',feature_names)
FileUtility.save_list('all_genomes_lens',[str(x) for x in [len(train),len(test),len(negative_class)]])
