# DiTaxa Installation
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

# Installation instruction using miniconda on Linux and MacOSx

We recommend to use a conda virtual environment for installation of DiTaxa.


## (1) Miniconda installation

The first step is to install the latest version of conda on your system.

### Linux
```
cd ~
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh
```

### MacOS
```
cd ~
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash ./Miniconda3-latest-MacOSX-x86_64.sh
```


## (2) General configurations of conda envinronment

Then you need to add the conda to your path. Please modify the following path if you changed in the installation

```
export PATH="~/miniconda3/bin:$PATH"
```

Then you need to add conda channels:

```
conda config --add channels conda-forge
conda config --add channels bioconda
```


## (3) Cloning the DiTaxa repository and install NCBI BLAST

Either download the zip file of DiTaxa from the <a href='https://github.com/ehsanasgari/DiTaxa/archive/master.zip'>github repository </a> or use the following command:

```
git clone git@github.com:ehsanasgari/DiTaxa.git
cd DiTaxa/
```


### NCBI BLAST

Run build.sh to download the latest blast alignment tool for your system (linux or MacOS). This bash script will automatically download it and move it in the default
directory that DiTaxa uses (DiTaxa/ncbi-blast/bin/). If you already have this tool in your system you can skip this part and just provide its path to the program later.

```
bash ./build.sh
```


## (4) Installation of dependencies in the virtual environment

The next step would be installation of the dependencies:

### Linux
```
conda create --name DiTaxa --file installations/env_linux.txt
```

### MacOS
```
conda create --name DiTaxa --file installations/env_macosx.txt

```

### Linux and MacOS


Then you need to activate the DiTaxa virtual environment:

```
source activate DiTaxa
```

Install the <a href='https://github.com/google/sentencepiece/tree/master/python
'>Google SentencePiece </a> via pip:

```
pip install sentencepiece
```


## (5) Testing the installation

For a quick test you may run the following command:

```
python ditaxa.py -h
```

In order to run the working example you can run:

```
bash ./run_test.sh
```

or its content, i.e.:

```
python3 ditaxa.py --indir dataset/periodontal/ --fast2label dataset/periodontal/mapping.txt --ext fastq --outdir results_dental/ --dbname periodontal --cores 4 --phenomap diseased:1,healthy:0 --phenoname DvsH --override 0 --heatmap periodontal:healthy

```

This should generate the output files (the intermediate representations, the biomarkers excel file, the heatmap, and the phylogenetic tree, etc.) step by step in the results_dental/ directory.
For more details on the input and the outputs please see <a href='https://github.com/ehsanasgari/DiTaxa'>here</a>.
