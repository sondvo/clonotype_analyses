# clonotype_analyses
**Notice**: *A python package is currently in active development at branch [api](https://github.com/sondvo/clonotype_analyses/tree/api)*

## Overview:
- This is a package for the analysis of single-cell TCR/BCR data, incorperating with single-cell RNA clinical metadata
- Functions include:
	- preprocessing
	- ratio of normal and ambigous clonotypes
	- shannon entropy
	- clonotype expansion
	- cdr3 length
	- clonotype diversity on embedding
	- clonotype tracing
- Every analysis will have charts from both [matplotlib](./TCR_analysis.ipynb) and [plotly](./html_plot.py)


### Example dataset: GSE185381
- Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185381
- Used datasets: xxx_vdj_t_**filtered**_contig_annotations.csv.gz
- Finish plotly and dash
- Viewing results:
`python3 plotly_html.py`

![alt text](plotly_html.gif)