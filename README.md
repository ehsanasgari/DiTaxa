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


<strong>Contact/Developer</strong>: Ehsaneddin Asgari (<span style="color: #0000ff;">asgari [at] berkeley [dot] edu</span>)
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

Please cite the <a style="color: #800000;" href="https://www.biorxiv.org/content/early/2018/01/31/255018">bioarXiv</a> version  <a href="https://www.biorxiv.org/highwire/citation/78275/bibtext"><img class="alignnone wp-image-142" src="http://llp.berkeley.edu/wp-content/uploads/2018/01/bibtex-icon.png" alt="" width="44" height="44" /></a> <a href="https://www.biorxiv.org/highwire/citation/78275/mendeley"><img class="alignnone wp-image-143" src="http://llp.berkeley.edu/wp-content/uploads/2018/01/Apps-Mendeley-icon-150x150.png" alt="" width="47" height="41" /></a>

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

DiTaxa is implemented in Python3.x and uses Scikit-Learn and Keras frameworks for machine learning. To install the dependencies use the following command:
```
pip install -r requirements.txt
```


<h2>Local ezCloud blast and GraPhlAn setup</h2>

The local blast files are already compiled and exist in directory "db".
Please make sure that NCBI BLASTN is already installed in your system.
If blastn is not added to the path, you can alternatively add the path to
"db/blastn_config.txt" file.

<h1> User Manual </h1>

```
python3 ditaxa.py --indir address_of_samples --ext extension_of_the_files --outdir output_directory --dbname database_name --cores 20 --fast2label mapping_file_from_name_to_phenotype --phenomap mapping_labels_to_binary_1_or_0_phenotype
```

Using the above mentioned command all the steps will be done sequentially and output will be organized in subdirectories.

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

<h1> Running Example </h1>

An example of periodontal disease dataset (Jorth et al, 2015) is provided. You may run the following command after installation:

```
python3 ditaxa.py --indir dataset/periodontal/
 --fast2label dataset/periodontal/mapping.txt
 --ext fastq
 --outdir ../results_dental/
 --dbname periodontal
 --cores 20
 --phenomap diseased:1,healthy:0
 --phenoname DvsH
 --override 1

```

You can see that "dataset/periodontal/" contains fastq files for each sample, "dataset/periodontal/mapping.txt" provides a mapping from fastq files to their labels.
"--phenomap diseased:1,healthy:0" determines which labels to be considered as positive class and which as negative class. "phenoname", i.e. DvsH defines the name of this
labeling scheme (disease versus healthy).

<h1> Output example </h1>




<h2>Bootstrapping for sample size selection</h2>
<img src="https://user-images.githubusercontent.com/8551117/40692939-f8b2785c-63b4-11e8-9194-c944775bbdf6.png">

<h2>Phenotype classification</h2>
Supporting Random Forest, SVM, Cross-entropy classifier, Neural netwrok

<h2>Biomarker detection</h2>

<img src="https://user-images.githubusercontent.com/8551117/40692994-568c1960-63b5-11e8-999d-ee69917943f5.png"/>


<h2>Heatmap creation</h2>

<img src="https://user-images.githubusercontent.com/8551117/40692895-b01a2a68-63b4-11e8-95d2-fc3727471bca.png"/>



