# Description

This GitHub repository contains the scripts and pipelines used to generate bioinformatic data related to the analysis of *Drosophila* homologs of CNV and neurodevelopmental genes in the fly wing.

There are three directories in this repository: 


1. K-means clustering of adult wing qualitative phenotype data.
 * This directory contains an R script to perform k-means clustering of wing phenotypic data, and a text file containing the input phenotypic scores for each RNAi line.

2. Tissue-specific network analysis of CNV and neurodevelopmental genes.
 * This directory contains a Python script (`nearest_neighbor_weighted_allgenes.py`) for analyzing the connectivity of genes within a human tissue-specific interaction networks. The script takes an individual gene (gene name and Entrez ID) and tissue type as standard input, and outputs the shortest distance to every other gene in the network, as well as the genes located in the shortest path between the two target genes. For example, to find the connectivity of NCBP2 (Entrez ID 22916) to all other genes in the genome in a brain-specific network, one would run:
 `python nearest_neighbor_weighted_allgenes.py NCBP2_22916 brain`
 * The text file within the directory (`genes.protein-coding.txt`) contains a set of all protein-coding genes and Entrez IDs, used in conjunction with the script.
 * The network files used in this analysis (i.e. `brain_top_min_0.2.txt`) are described in [Greene *et al, Nat. Genet.* 2015](https://www.ncbi.nlm.nih.gov/pubmed/25915600), and can be obtained from the open-access GIANT database within [HumanBase](https://hb.flatironinstitute.org/download); we extracted all edges from each "Top Edges" network file from HumanBase with edge weight >0.2 to create a sub-network of highly connected genes for our analysis. 

3. Phenotype ontology analysis for carriers of pathogenic CNVs.
* This directory contains a Python script (`convert_hpo.py`) that takes a file listing human phenotypes (in our case, phenotypes of pathogenic CNV carriers from the Decipher database) and summarizes each phenotype
into its top-level Human Phenotype Ontology term. The script can be run using `python convert_hpo.py file_name`.
* The input file is a tab-delimited text file with a pathogenic CNV name in the 1st column, sample ID in the 2nd column, and comma-separated list of phenotypes in the 10th column. These files were curated from data available in the open-access [Decipher database](https://decipher.sanger.ac.uk/).
* The script also requires the most recent Human Phenotype Ontology .obo file, which can be downloaded from the [Human Phenotype Ontology website](http://purl.obolibrary.org/obo/hp.obo).
* The output file is a tab-separated text file listing the CNV name, sample ID, input phenotype list, comma-separated lists of HPO IDs for each input phenotype, and tab-separated lists of top-level HPO terms.
* If a phenotype is not present within HPO, the script will print the phenotype to standard output in an error message. The user can then manually modify the phenotypes to a similar term within HPO and re-run the script. (The script recognizes and skips phenotype modifiers, such as "moderate", "severe", or "family predisposition").

R scripts can be run using any R version (scripts were generated using R v.3.6.1.). Python scripts for network analysis can be run in Python2 (scripts were generated using Python v.2.7.16) and require the [NetworkX package v.2.4](https://networkx.github.io/) and the [Orange3 Bioinformatics package](https://orange-bioinformatics.readthedocs.io/). 

# Citation
Yusuff T, Jensen M, Yennawar S, Pizzo L, Karthikeyan S, Gould DJ, Sarker A, Matsui Y, Iyer J, Lai ZC, Girirajan S. *Drosophila* models of pathogenic copy-number variant genes show global and non-neuronal defects during development. 

# Copyright/License
The code in this repository is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the code.  If not, see <https://www.gnu.org/licenses/>.

# Contact
For questions or comments, please contact Matthew Jensen (mpj5142@psu.edu) or Santhosh Girirajan (sxg47@psu.edu).