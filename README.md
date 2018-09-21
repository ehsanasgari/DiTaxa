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
 Identifying distinctive taxa for microbiome-related diseases is considered key to the establishment of diagnosis and therapy options in precision medicine and imposes high demands on the accuracy of microbiome analysis techniques. We propose an alignment- and reference- free subsequence based 16S rRNA data analysis,
 as a new paradigm for microbiome phenotype and biomarker detection. Our method, called DiTaxa, substitutes standard OTU-clustering by segmenting 16S rRNA reads into the most frequent variable-length subsequences. We compared the performance of DiTaxa to the state-of-the-art methods in phenotype and biomarker detection,
 using human-associated 16S rRNA samples for periodontal disease, rheumatoid arthritis, and inflammatory bowel diseases, as well as a synthetic benchmark dataset. DiTaxa performed competitively to the k-mer based state-of-the-art approach in phenotype prediction while outperforming the OTU-based state-of-the-art approach
 in finding biomarkers in both resolution and coverage evaluated over known links from literature and synthetic benchmark datasets.
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

For the detailed installation using conda virtual environment and testing the working example please refer to the <a href='https://github.com/ehsanasgari/DiTaxa/tree/master/installations'> installation guideline </a>.

<h1> Working Example </h1>

An example of periodontal disease dataset (Jorth et al, 2015) is provided in the repo. In order to see how DiTaxa runs, you may run the following command after <a href='https://github.com/ehsanasgari/DiTaxa/tree/master/installations'>installation</a>.
<br/>
You need to replace BLASTN_PATH with the latest version of blast for your operating system, which you can get from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/):

```
python3 ditaxa.py --indir dataset/periodontal/
 --fast2label dataset/periodontal/mapping.txt
 --ext fastq
 --outdir results_dental/
 --dbname periodontal
 --cores 20
 --phenomap diseased:1,healthy:0
 --heatmap PeriodontalSamples:HealthySamples
 --phenoname DvsH
 --override 1
 (optional)--blastn BLASTN_PATH
```

Alternatively you can run:
```
bash ./run_test.sh
```

<h2> Example Dataset and parameter explanation</h2>
You may use this example to prepare your input files:

The "indir": e.g. «dataset/periodontal/» contains fastq files for each 16S rRNA samples.<hr/>
The "fast2label"" e.g. «dataset/periodontal/mapping.txt» provides a file containing mapping from fastq files to their labels in a tabular format:
```
d1.fastq    diseased
d2.fastq    diseased
d3.fastq    diseased
d4.fastq    diseased
d5.fastq    diseased
d6.fastq    diseased
d7.fastq    diseased
d8.fastq    diseased
d9.fastq    diseased
d10.fastq    diseased
h1.fastq    healthy
h2.fastq    healthy
h3.fastq    healthy
h4.fastq    healthy
h5.fastq    healthy
h6.fastq    healthy
h7.fastq    healthy
h8.fastq    healthy
h9.fastq    healthy
h10.fastq    healthy
```
<hr/>
The "phenomap", e.g. «diseased:1,healthy:0» determining which labels to be considered as positive class and which as negative class as a string with no space in the following format:

```
diseased:1,healthy:0
```

<hr/>
The "override", 1 will override already existing files in the directory.<hr/>
The "heatmap" e.g. «PeriodontalSamples:HealthySamples» determines the names for plotting positive and negative pheotypes on the heatmap.<hr/>
The "blastn", (optional: only if you don't run build.sh you need to specify this) path to the "bin" directory of blast existing on your system.
<br/>

After running this command the output files will be generated in 'results_dental' as described bellow. The example output files are provided in the './output_example/'
directory.

<h2> Output example </h2>
The automatically generated output of the example is as follows, you may also see automatically generated files in `output_example`:

![ditaxaout](https://user-images.githubusercontent.com/8551117/42161908-1c2d9d72-7dfd-11e8-86eb-68af5055c5ab.png)

<h3>Biomarker detection</h3>
DiTaxa provides a taxonomic tree for significant discriminative biomarkers, where identified taxa to the positive and negative class are colored according to their phenotype (red for positive class and blue for negative class). The DiTaxa implementation for taxonomic tree generation uses a Phylophlan-based backend.

<img src="https://user-images.githubusercontent.com/8551117/40692994-568c1960-63b5-11e8-999d-ee69917943f5.png"/>


<h3>Heatmap of Biomarkers</h3>

DiTaxa provides a heatmap of top biomarkers occurrences in samples, where the rows denote markers and the columns are samples is generated. Such a heatmap allows biologists to obtain a detailed overview of markers' occurrences across samples. The heatmap shows number of distinctive sequences hit by each biomarker in different samples and stars in the heatmap denote hitting unique sequences, which cannot be analyzed by OTU clustering approaches.

<img src="https://user-images.githubusercontent.com/8551117/40692895-b01a2a68-63b4-11e8-95d2-fc3727471bca.png"/>


<h3>Excel file of Biomarkers</h2>
In addition, DiTaxa provides a detailed excel file of biomarker sequnces and their taxonomy annotations along with their p-values.

<h3>t-SNE visualization</h3>
T-sne visualization of data using all NPEs and selected markers will be also generated by default.




<h1> Detailed User Manual </h1>

After installation using the <a href='https://github.com/ehsanasgari/DiTaxa/tree/master/installations'> installation guideline </a>.you may use DiTaxa The parameteres for running DiTaxa are as follows:

```
python3 ditaxa.py --indir address_of_samples --ext extension_of_the_files --outdir output_directory --dbname database_name --cores 20 --fast2label mapping_file_from_name_to_phenotype --phenomap mapping_labels_to_binary_1_or_0_phenotype
--blastn /mounts/data/proj/asgari/dissertation/deepbio/taxonomy/ncbi-blast-2.5.0+/bin/

```

Using the above mentioned command all the steps will be done sequentially and output will be organized in subdirectories.

<h3> Main parameters for biomarker detection/analysis </h3>

--indir: The input directory containing all fasta or fastq files. (e.g.: datasets/periodontal/)<br/>
--ext: Sequence file extensions (fasta or fastq) (e.g.: fastq)<br/>
--outdir: The output directory (e.g.: /mounts/data/ditaxa/results/test_dental_out/)<br/>
--cores: Number of cores (e.g.: 40)<br/>
--fast2label: tabular mapping file between file names and the labels<br/>
--phenomap: mapping from label to binary phenotypes<br/>
--phenoname: name of the phenotype mapping, if not given the labels and their value will be used for identification: label1@1#label2@1...#label3@0. Please note that a single project may have several phenotype mapping schemes (untreated diseased versus all or untreated versus healthy or etc.)<br/>
--override: 1 to override the existing files, 0 to only generate the missing files<br/>
--heatmap: generates occurrence heatmap for the top 100 markers (e.g:  positive_title:negative_title).<br/>
--excel: 1 or 0, the default is 1 to generate a detailed list of markers, their taxonomic assignment, and their p-values<br/>
--blastn:  If you have already run './build.sh' you do not need to specify this parameter and the script will download it and put the
<a href=' ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/'>NCBI BLASTN</a> /bin/ path in your system. Otherwise, if you already have this on your system you can specify it here.
<p>
You can also download blast+ from below and specify the path:<br/>
<b>Linux</b><br/>
http://ftp.ncbi.nlm.nih.gov/blast/executables/blast%2B/2.7.1/ncbi-blast-2.7.1%2B-x64-linux.tar.gz
<br/><b>MacOSx</b><br/>
http://ftp.ncbi.nlm.nih.gov/blast/executables/blast%2B/2.7.1/ncbi-blast-2.7.1%2B-x64-macosx.tar.gz
</p>

<h3> Phenotype prediction </h3>

For phenotype classification functionality, evaluation a 10XFold cross-validation framework:
--classify: which predictive model to use: choices=[False: default, 'RF': random forest, 'SVM': support vector machines, 'DNN': deep multi-layer perceptron, 'LR': logistic regression] <br/>

<b>Deep neural network parameters</b>

Although a full script is provided, in order to simplify the core installation of DiTaxa for biomarker detection/analysis we have commented the deep neural network classifier and its dependencies. In case you are interested in using neural network prediction of the phenotype you only need to install some further dependencies (keras/tensorflow) and uncomment "import DNN" in main/DiTaxa.py.

--arch: The comma separated definition of neural network layers connected to eahc other, you do not need to specify the input and output layers, values between 0 and 1 will be considered as dropouts, e.g., 1024,0.2,512'<br/>
--batchsize<br/>
--gpu_id: which GPU to use<br/>
--epochs: Number of epochs<br/>

<h3>Bootstrapping for sample size selection</h3>
We use bootstrapping to investigate sufficiency and consistency of NPE representation, when only a small portion of the sequences are used. This has two important implications, first, sub-sampling reduces the preprocessing run-time, second, it shows that even a shallow 16S rRNA sequencing is enough for the phenotype prediction. We use a resampling framework to find a proper sampling size. The DiTaxa implementation uses a defualt parameter setting based on bootstrapping on several datasets. The bootstrapping library is located at "DiTaxa/bootstrapping/bootstrapping.py" if further investigation is needed.

<img src="https://user-images.githubusercontent.com/8551117/40692939-f8b2785c-63b4-11e8-9194-c944775bbdf6.png">



