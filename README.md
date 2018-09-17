# DiTaxa
<table style="height: 48px; width: 812px;">
<table style="width: 802px;">
<tbody>
<tr>
<td style="width: 450px;" colspan="2"><span style="font-size: 40pt; font-family: helvetica,arial,sans-serif;"><span style="color: #0000ff;"><strong>Nucleotide-pair encoding of 16S rRNA sequences for host phenotype and biomarker detection</strong></span></span></td>
</tr>
<tr>
<td style="width: 450px;background-color: white;" colspan="2"><img src="https://user-images.githubusercontent.com/8551117/40691993-7b6e0014-63af-11e8-9289-c3b842aff9f6.png"/></td>
</tr>
</tbody>
</table>


Asgari E., Münch P. C., Lesker T. R., McHardy A. C. &#9733;, and Mofrad M. R. K. &#9733;,<br/>
<b>Nucleotide-pair encoding of 16S rRNA sequences for host phenotype and biomarker detection.</b><br/>
bioRxiv (2018): 334722.


<strong>Developer</strong>: Ehsaneddin Asgari (<span style="color: #0000ff;">asgari [at] berkeley [dot] edu</span>)
<br/>
Please feel free to report any technical issue by sending an email or reporting an issue here.
<br/>
<strong>Project page:</strong> <a href="http://llp.berkeley.edu/ditaxa">http://llp.berkeley.edu/ditaxa</a>
<br/>
<strong>PIs</strong>: Prof. Alice McHardy* and Prof. Mohammad Mofrad*
<br/>
 <img class="alignnone wp-image-125" src="http://llp.berkeley.edu/wp-content/uploads/2018/01/Logo_HZI_2x1-300x129.png" alt="" width="202" height="87" /> <img class="alignnone size-full wp-image-9" src="http://llp.berkeley.edu/wp-content/uploads/2018/01/logo_ucb2-e1516672757796.png" alt="" width="100" height="97" />       

<hr />

<span style="font-family: helvetica,arial,sans-serif; font-size: 24pt;"><strong>Summary</strong></span>

&nbsp;
Identifying combinations of taxa distinctive for microbiome-associated diseases is considered key to the establishment of diagnosis and therapy options in precision medicine and imposes high demands on accuracy of microbiome analysis techniques. We propose subsequence based 16S rRNA data analysis, as a new paradigm for microbiome phenotype classification and biomarker detection. This method and software called DiTaxa substitutes standard OTU-clustering or sequence-level analysis by segmenting 16S rRNA reads into the most frequent variable-length subsequences. These subsequences are then used as data representation for downstream phenotype prediction, biomarker detection and taxonomic analysis. Our proposed sequence segmentation called nucleotide-pair encoding (NPE) is an unsupervised data-driven segmentation inspired by Byte-pair encoding, a data compression algorithm. The identified subsequences represent commonly occurring sequence portions, which we found to be distinctive for taxa at varying evolutionary distances and highly informative for predicting host phenotypes.
We compared the performance of DiTaxa to the state-of-the-art methods in disease phenotype prediction and biomarker detection, using human-associated 16S rRNA samples for periodontal disease, rheumatoid arthritis and inflammatory bowel diseases, as well as a synthetic benchmark dataset. DiTaxa identified 17 out of 29 taxa with confirmed links to periodontitis (recall=0.59), relative to 3 out of 29 taxa (recall=0.10) by the state-of-the-art method. On synthetic benchmark data, DiTaxa obtained full precision and recall in biomarker detection, compared to 0.91 and 0.90, respectively. In addition, machine-learning classifiers trained to predict host disease phenotypes based on the NPE representation performed competitively to the state-of-the art using OTUs or k-mers. For the rheumatoid arthritis dataset, DiTaxa substantially outperformed OTU features with a macro-F1 score of 0.76 compared to 0.65. Due to the alignment- and reference free nature, DiTaxa can efficiently run on large datasets. The full analysis of a large 16S rRNA dataset of 1359 samples required ~1.5 hours on 20 cores, while the standard pipeline needed ~6.5 hours in the same setting.
&nbsp;</td>
</tr>
</tbody>

</table>

Please cite the <a style="color: #800000;" href="https://www.biorxiv.org/content/early/2018/05/30/334722">bioarXiv</a> version  <a href="https://www.biorxiv.org/highwire/citation/101986/bibtext"><img class="alignnone wp-image-142" src="http://llp.berkeley.edu/wp-content/uploads/2018/01/bibtex-icon.png" alt="" width="44" height="44" /></a> <a href="https://www.biorxiv.org/highwire/citation/101986/mendeley"><img class="alignnone wp-image-143" src="http://llp.berkeley.edu/wp-content/uploads/2018/01/Apps-Mendeley-icon-150x150.png" alt="" width="47" height="41" /></a>

```
@article{asgari2018nucleotide,
  title={Nucleotide-pair encoding of 16S rRNA sequences for host phenotype and biomarker detection},
  author={Asgari, Ehsaneddin and M{\"u}nch, Philipp C and Lesker, Till R and McHardy, Alice Carolyn and Mofrad, Mohammad RK},
  journal={bioRxiv},
  pages={334722},
  year={2018},
  publisher={Cold Spring Harbor Laboratory}
}

```



<h1>Installation</h1>

DiTaxa implementation in Python3.6.x is provided (we recommend 3.6.5). You may use <a href='http://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/'>virtual environment</a>.

To install the dependencies use the following command:
```
pip3 install --upgrade pip
pip3 install -r requirements.txt
```
or:
```
conda install --yes --file requirements.txt
```

DiTaxa for segmentation of the sequences by default uses Google SentencePiece backend. For non-Linux environment, including macOS, Windows, and Linux (arm) please make sure that you  have installed SentencePiece C++ library in advance as instructed here:

```
https://github.com/google/sentencepiece/tree/master/python
```

<h1> User Manual </h1>

```
python3 ditaxa.py --indir address_of_samples --ext extension_of_the_files --outdir output_directory --dbname database_name --cores 20 --fast2label mapping_file_from_name_to_phenotype --phenomap mapping_labels_to_binary_1_or_0_phenotype
--blastn /mounts/data/proj/asgari/dissertation/deepbio/taxonomy/ncbi-blast-2.5.0+/bin/

```

Using the above mentioned command all the steps will be done sequentially and output will be organized in subdirectories.

<h3> Main parameters and biomarker detection </h3>

--indir: The input directory containing all fasta or fastq files. (e.g.: datasets/periodontal/)<br/>
--ext: Sequence file extensions (fasta or fastq) (e.g.: fastq)<br/>
--outdir: The output directory (e.g.: /mounts/data/ditaxa/results/test_dental_out/)<br/>
--cores: Number of cores (e.g.: 40)<br/>
--fast2label: tabular mapping file between file names and the labels<br/>
--phenomap: mapping from label to binary phenotypes<br/>
--phenoname: name of the phenotype mapping, if not given the labels and their value will be used for identification: label1@1#label2@1...#label3@0. Please note that a single project may have several phenotype mapping schemes (untreated diseased versus all or untreated versus healthy or etc.)<br/>
--override: 1 to override the existing files, 0 to only generate the missing files<br/>
--heatmap: generates occurrence heatmap for top 100 markers (e.g:  positive_title:negative_title).<br/>
--excel: 1 or 0, the default is 1 to generate a detailed list of markers, their taxonomic assignment, and their p-values<br/>
--blastn:  NCBI BLASTN path in your system, you can get the latest version from here: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/


<h3> Phenotype prediction </h3>
For phenotype classification functionality, evaluation a 10XFold cross-validation framework:
--classify: which predictive model to use: choices=[False: default, 'RF': random forest, 'SVM': support vector machines, 'DNN': deep multi-layer perceptron, 'LR': logistic regression] <br/>

<b>Deep neural network parameters</b>
--arch: The comma separated definition of neural network layers connected to eahc other, you do not need to specify the input and output layers, values between 0 and 1 will be considered as dropouts, e.g., 1024,0.2,512'<br/>
--batchsize<br/>
--gpu_id: which GPU to use<br/>
--epochs: Number of epochs<br/>


<h1> Working Example </h1>

An example of periodontal disease dataset (Jorth et al, 2015) is provided (a relatively small dataset). In order to see how DiTaxa runs, you may run the following command after installation. You need to replace BLASTN_PATH with the latest version of blast for your operating system, which you can get from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/):

```
python3 ditaxa.py --indir dataset/periodontal/
 --fast2label dataset/periodontal/mapping.txt
 --ext fastq
 --outdir results_dental/
 --dbname periodontal
 --cores 20
 --phenomap diseased:1,healthy:0
 --phenoname DvsH
 --override 1
 --blastn BLASTN_PATH
```

<h2> Example Dataset </h2>
You may use this example to prepare your input files:

The "indir": e.g. «dataset/periodontal/» contains fastq files for each sample.<br/>
The "fast2label"" e.g. «dataset/periodontal/mapping.txt» provides a mapping from fastq files to their labels.<br/>
The "phenomap", e.g. «diseased:1,healthy:0» determining which labels to be considered as positive class and which as negative class.<br/>
The "override", 1 will override already existing files in the directory.<br/>
The "blastn", path to the "bin" directory of blast existing on your system.
You can get the latest version from here: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
<br/>
After running this command the output files will be generated in 'results_dental' as described bellow.

<h2> Output example </h2>
The automatically generated output of the example is as follows:

![ditaxaout](https://user-images.githubusercontent.com/8551117/42161908-1c2d9d72-7dfd-11e8-86eb-68af5055c5ab.png)


<h2>Biomarker detection</h2>
DiTaxa provides a taxonomic tree for significant discriminative biomarkers, where identified taxa to the positive and negative class are colored according to their phenotype (red for positive class and blue for negative class). The DiTaxa implementation for taxonomic tree generation uses a Phylophlan-based backend.

<img src="https://user-images.githubusercontent.com/8551117/40692994-568c1960-63b5-11e8-999d-ee69917943f5.png"/>


<h2>Heatmap of Biomarkers</h2>

DiTaxa provides a heatmap of top biomarkers occurrences in samples, where the rows denote markers and the columns are samples is generated. Such a heatmap allows biologists to obtain a detailed overview of markers' occurrences across samples. The heatmap shows number of distinctive sequences hit by each biomarker in different samples and stars in the heatmap denote hitting unique sequences, which cannot be analyzed by OTU clustering approaches. 

<img src="https://user-images.githubusercontent.com/8551117/40692895-b01a2a68-63b4-11e8-95d2-fc3727471bca.png"/>

In addition, DiTaxa provides a detailed excel file of biomarker sequnces and their taxonomy annotations along with their p-values. T-sne visualization of data using all NPEs and selected markers will be also generated by default.


<h2>Bootstrapping for sample size selection</h2>
We use bootstrapping to investigate sufficiency and consistency of NPE representation, when only a small portion of the sequences are used. This has two important implications, first, sub-sampling reduces the preprocessing run-time, second, it shows that even a shallow 16S rRNA sequencing is enough for the phenotype prediction. We use a resampling framework to find a proper sampling size. The DiTaxa implementation uses a defualt parameter setting based on bootstrapping on several datasets. The bootstrapping library is located at "DiTaxa/bootstrapping/bootstrapping.py" if further investigation is needed.

<img src="https://user-images.githubusercontent.com/8551117/40692939-f8b2785c-63b4-11e8-9194-c944775bbdf6.png">



