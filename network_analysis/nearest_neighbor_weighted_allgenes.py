#This script takes an individual gene as standard input (gene name and Entrez ID), and outputs the shortest path length to every other gene in the network, 
#as well as the genes located in the shortest path between the two target genes. 
#For example, to find the connectivity of NCBP2 (Entrez ID 22916) to all other genes in the genome in the brain-specific network, one would run:
#python nearest_neighbor_weighted_allgenes.py NCBP2_22916 brain

#Import modules
import csv
import sys
import networkx
import itertools

#Open network file (format: Entrez_ID_1, Entrez_ID_2, Weighted_connectivity) and list of genes in network (Entrez IDs in column 1, gene names in column 3)
#The name of teh network to open comes from the command line
tissue_type=sys.argv[2]
network_file=open(tissue_type+"_top_min_0.2.txt",'rU')
network_lines=network_file.readlines()
gene_file=open("genes.protein-coding.txt",'rU')
gene_lines=gene_file.readlines()

#Populate list of protein-coding genes in the network, and dictionary for Entrez-gene symbol conversion
coding_genes_list=[]
entrez_dict={}
for gene_line in gene_lines:
	gene_line=gene_line.strip().split("\t")
	gene_ID=gene_line[0]
	gene_symbol=gene_line[2]
	coding_genes_list.append(gene_ID)
	entrez_dict[gene_ID]=gene_symbol
coding_genes=set(coding_genes_list)

#Get starting gene(s) from system input, and create list of target genes (i.e. every other gene in network) for each start gene
gene_targets={}
gene_input=sys.argv[1].split("_")
gene_targets[gene_input[0]]=[gene_input[1]]+coding_genes_list

#Populate network data for each gene with NetworkX
network_graph=networkx.Graph()
for network_line in network_lines:
	network_line=network_line.strip().replace('\x00','').split("\t")
	node1=network_line[0]
	node2=network_line[1]
	weight=1/float(network_line[2]) #Take inverse of weight for edge lengths

	#Add edge to NetworkX network if both genes are in protein-coding list	
	if node1 in coding_genes:
		if node2 in coding_genes:
			network_graph.add_edge(node1,node2,weight=weight)

#Iterate through list of start genes to test
for gene_name in gene_targets:

	#Initialize output file
    outwriter = csv.writer(file(gene_name+"_"+tissue_type+"_nearest_neighbor_weighted.txt",'w'),delimiter='\t')
    header=["Start gene","End gene","Shortest distance","Number of connector genes","Connector genes"]
	outwriter.writerow(header)

	#Get list of target genes for start gene (i.e. all other genes in network)
	gene_target_list=gene_targets[gene_name]
	gene_start_target=gene_target_list[0]
	gene_end_targets=gene_target_list[1:]

	#Iterate through all target genes
	for gene_end_target in gene_end_targets:
		print gene_start_target+"\t"+gene_end_target

		#Use NetworkX to generate the shortest paths between start and target genes, using Dijkstra's algorithm based on the weighted edges. 
		#Save the list of genes in the shortest path, the length of the shortest path (based on edge weights), and the number of connector genes in each path.
		try:
			shortest_path=list(networkx.dijkstra_path(network_graph,gene_start_target,gene_end_target))
			shortest_path_length=networkx.dijkstra_path_length(network_graph,source=gene_start_target,target=gene_end_target)
			num_conn_genes=len(shortest_path)-2
		except(networkx.exception.NetworkXNoPath): #Keep all fields empty if no path exists between the two genes.
			shortest_path=["No paths"]
			shortest_path_length=0
			num_conn_genes=0

		#Convert Entrez IDs of gene targets and genes in the shortest path to gene symbols
		gene_start_target_new=entrez_dict[gene_start_target] #Start gene
		try:
			gene_end_target_new=entrez_dict[gene_end_target] #Target gene
		except(KeyError):
			gene_end_target_new=gene_end_target #Skip if not in dictionary
		connector_genes_new=[]
		if shortest_path[0]!="No paths":
			for connector_gene in shortest_path: #Genes in shortest path
				try:
					connector_gene_new=entrez_dict[connector_gene]
				except(KeyError):
					connector_gene_new=connector_gene #Skip if not in dictionary
				connector_genes_new.append(connector_gene_new)
		connector_genes_new=connector_genes_new[1:-1] #Remove CNV and target gene from list of genes in shortest path

		#Output to file
		outrow=[gene_start_target_new,gene_end_target_new,shortest_path_length,num_conn_genes]+connector_genes_new
		outwriter.writerow(outrow)