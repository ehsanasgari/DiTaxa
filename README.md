# DiTaxa
<table style="height: 48px; width: 812px;">
<table style="width: 802px;">
<tbody>
<tr>
<td style="width: 450px;" colspan="2"><span style="font-size: 40pt; font-family: helvetica,arial,sans-serif;"><span style="color: #0000ff;"><strong>Nucleotide-pair encoding of 16S rRNA sequences for host phenotype and biomarker detection</strong></span></span></td>
</tr>
<tr>
<td style="width: 450px;background-color: white;" colspan="2"><img src="https://user-images.githubusercontent.com/8551117/40687526-355d0516-639b-11e8-894c-d843a1699433.png"/></td>
</tr>
</tbody>
</table>


Asgari E., Muench P., Lesker T.R., McHardy A.C. and Mofrad M.R.K., Nucleotide-pair encoding of 16S rRNA sequences for host phenotype and biomarker detection. bioRxiv, 2018. Available at: ======


 
The datasets </strong> are also available for download <a href='http://llp.berkeley.edu/ditaxa'><img class="alignnone wp-image-36" src="http://llp.berkeley.edu/wp-content/uploads/2018/01/zip.png" alt="" width="33" height="33" /></a>.

<strong>Contact</strong>: Ehsaneddin Asgari (<span style="color: #0000ff;">asgari [at] berkeley [dot] edu</span>)
<br/>
<strong>Project page:</strong>: <a href="http://llp.berkeley.edu/ditaxa">http://llp.berkeley.edu/ditaxa</a>
<hr />

<span style="font-family: helvetica,arial,sans-serif; font-size: 24pt;"><strong>Summary</strong></span>

&nbsp;

We propose subsequence based 16S rRNA data processing, as a new paradigm for sequence phenotype classification and biomarker detection. This method and software called DiTaxa substitutes standard OTU-clustering or sequence-level analysis by segmenting 16S rRNA reads into the most frequent variable-length subsequences. These subsequences are then used as data representation for downstream phenotype prediction, biomarker detection and taxonomic analysis. Our proposed sequence segmentation called nucleotide-pair encoding (NPE) is an unsupervised data-driven segmentation inspired by Byte-pair encoding, a data compression algorithm. The identified subsequences represent commonly occurring sequence portions, which we found to be distinctive for taxa at varying evolutionary distances and highly informative for predicting host phenotypes.
We compared the performance of DiTaxa to the state-of-the-art methods in disease phenotype prediction and biomarker detection, using human-associated 16S rRNA samples for periodontal disease, rheumatoid arthritis and inflammatory bowel diseases, as well as a synthetic benchmark dataset. DiTaxa identified 13 out of 21 taxa with confirmed links to periodontitis (recall=0.62), relative to 3 out of 21 taxa (recall=0.14) by the state-of-the-art method. On synthetic benchmark data, DiTaxa obtained full precision and recall in biomarker detection, compared to 0.91 and 0.90, respectively. In addition, machine-learning classifiers trained to predict host disease phenotypes based on the NPE representation performed competitively to the state-of-the art using OTUs or k-mers. For the rheumatoid arthritis dataset, DiTaxa substantially outperformed OTU features with a macro-F1 score of 0.76 compared to 0.65. Due to the alignment- and reference free nature, DiTaxa can efficiently run on large datasets. The full analysis of a large 16S rRNA dataset of 1359 samples required ~1.5 hours on 20 cores, while the standard pipeline needed ~6.5 hours in the same setting.

&nbsp;</td>
</tr>
</tbody>

</table>


<h1>Installation</h1>

MicroPheno is implemented in Python3.x and uses ScikitLearn and Keras frameworks for machine learning. To install the dependencies use the following command:
```
pip install -r requirements.txt
```

Please cite the <a style="color: #800000;" href="https://www.biorxiv.org/content/early/2018/01/31/255018">bioarXiv</a> version  <a href="https://www.biorxiv.org/highwire/citation/78275/bibtext"><img class="alignnone wp-image-142" src="http://llp.berkeley.edu/wp-content/uploads/2018/01/bibtex-icon.png" alt="" width="44" height="44" /></a> <a href="https://www.biorxiv.org/highwire/citation/78275/mendeley"><img class="alignnone wp-image-143" src="http://llp.berkeley.edu/wp-content/uploads/2018/01/Apps-Mendeley-icon-150x150.png" alt="" width="47" height="41" /></a>

```
@article {Asgari255018,
	author = {Asgari, Ehsaneddin and Garakani, Kiavash and McHardy, Alice Carolyn and Mofrad, Mohammad R.K.},
	title = {MicroPheno: Predicting environments and host phenotypes from 16S rRNA gene sequencing using a k-mer based representation of shallow sub-samples},
	year = {2018},
	doi = {10.1101/255018},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2018/01/31/255018},
	eprint = {https://www.biorxiv.org/content/early/2018/01/31/255018.full.pdf},
	journal = {bioRxiv}
}

```

<h1> User Manual </h1>
MicroPheno can be used either via the templates provided in the <a href="https://github.com/ehsanasgari/MicroPheno/tree/master/notebooks">ipython notebooks</a> or the command-line interface.

<h2>Bootstrapping</h2>
An example of bootstrapping provided in the <a href="https://github.com/ehsanasgari/MicroPheno/blob/master/notebooks/1.Bootstrapping.ipynb">notebooks</a>.

<b>Command line use:</b> Argument to be used are the input/output directories, the sequence filetype, the k-mers and the sample size. Use argument '-h' to see the helpers.
```
python3 micropheno.py --bootstrap --indir /path/to/16srRNAsamples/ --out output_dir/ --filetype fastq --kvals 3,4,5,6 --nvals 10,100,200,500,1000 --name crohs
```
The output would be generating the following plot in the specified output directory. See the related <a href="https://github.com/ehsanasgari/MicroPheno/blob/master/notebooks/1.Bootstrapping.ipynb">notebook</a> for more details.
<img src="https://user-images.githubusercontent.com/8551117/35446008-af953ad6-02b3-11e8-9b33-06d1f4b429f3.png" alt="bootstrapping" />


<h2>Representation Creation</h2>
Two examples of representation creation are provided in the <a href="https://github.com/ehsanasgari/MicroPheno/blob/master/notebooks/2.%20k-mer%20Representation%20Creation%20with%20sub-sampling%20or%20without.ipynb">notebooks</a>, one with sampling from sequence files and the other for mapping the representative sequences.

<b>Command line use:</b> Argument to be used are the input/output directories, the sequence filetype, the k-mers and their sample size as well as number of cores to be used. Use argument '-h' to see the helpers.

```
python3 micropheno.py --genkmer --inaddr /path/to/16srRNAsamples/ --out output_dir/ --filetype fastq --cores 20 --KN 6:100,6:1000,2:100 --name test_crohn
```

<h2>Classification with Random Forest and SVM</h2>

The trained representation in the previous step in the input for classification.
See an example in the<a href="https://github.com/ehsanasgari/MicroPheno/blob/master/notebooks/3.%20Classification_classical_classifiers.ipynb"> notebooks</a>.

<b>Command line use:</b> Argument to be used are the X and Y, the classification algorithm (SVM, or RF), output directory as well as number of cores to be used. Use argument '-h' to see the helpers.

The following command will do tuning the parameters as well as evaluation within a 10xFold corss-validation scheme. Details on how to parse the results (scores, confusion matrix, best estimator, etc) are provided <a href="https://github.com/ehsanasgari/MicroPheno/blob/master/notebooks/3.%20Classification_classical_classifiers.ipynb"> here</a>.

```
python3 micropheno.py --train_predictor --model RF (or SVM) --x k-mer.npz --y labels_phenotypes.txt --cores 20 --name test_crohn  --out output_dir/
```

<h2>Classification with Deep Neural Network</h2>
We use the Multi-Layer-Perceptrons (MLP) Neural Network architecture with several hidden layers using Rectified Linear Unit (ReLU) as the nonlinear activation function. We use softmax activation function at the last layer to produce the probability vector that can be regarded as representing posterior probabilities (Goodfellow-et-al-2016). To avoid overfitting we perform early stopping and also use dropout at hidden layers (Srivastava2014). A schematic visualization of our Neural Networks is depicted in the Figure.

<img src="https://user-images.githubusercontent.com/8551117/35446216-4ec1eb7c-02b4-11e8-9421-043ec1f9ed96.png" alt="dnn" />

Our objective is minimizing the loss, i.e. cross entropy between output and the one-hot vector representation of the target class. The error (the distance between the output and the target) is used to update the network parameters via a Back-propagation algorithm using Adaptive Moment Estimation (Adam) as optimizer (Kingma2015).

You can see an example in the notebooks <a href="https://github.com/ehsanasgari/MicroPheno/blob/master/notebooks/4.%20Classification%20Deep%20Learning.ipynb">here</a>, showing how to see the learning curves and also getting the activation function of the neural network from the trained model.

<b>Command line use:</b> Argument to be used are the X and Y, the DNN flag, the neural network architecture (hidden-sizes and dropouts), batch size, number of epochs, output directory as well as the GPU id to be used. Use argument '-h' to see the helpers.

```
python3 micropheno.py --train_predictor --model DNN --arch 1024,0.2,512,0.1,128,64 --batchsize 10 --epochs  100 --x k-mer.npz --y labels_phenotypes.txt --name test_crohn  --out output_dir/
```


<h2>Visualization</h2>

An example of visualization using PCA, t-SNE, as well as t-SNE over the activation function of the last layer of the neural network is provided in <a href="https://github.com/ehsanasgari/MicroPheno/blob/master/notebooks/5.%20Visualization.ipynb">this notebook</a>.


![vis](https://user-images.githubusercontent.com/8551117/35447281-8f58b064-02b7-11e8-9a97-affe35573ba5.png)


