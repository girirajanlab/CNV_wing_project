#This Python script takes a list of human phenotypes (in our case, phenotypes of pathogenic CNV carriers from the Decipher database) and summarizes each phenotype
#into its top-level Human Phenotype Ontology term, using the Orange3 Bioinformatics library to handle ontology data structures. 
#The input file is a tab-delimited text file with a pathogenic CNV name in the 1st column, sample ID in the 2nd column, and comma-separated list of phenotypes in the 10th column. 
#The output file is a tab-separated text file listing the CNV name, sample ID, input phenotype list, comma-separated lists of HPO IDs for each input phenotype, and tab-separated lists of top-level HPO terms. 
#The script also requires the most recent Human Phenotype Ontology .obo file, names as "hp.obo" below.
#The script can be called as: "python convert_hpo.py file_name"

#Method for identifying top-level HPO term (input: HPO term name and existing list of top-level terms)
def getNextTerms(name,top_level):
	top_level_list=top_level
	term=ontology.term(name) #Look up term in ontology
	related_terms=term.related_objects()
	next_terms=[rel_term[1] for rel_term in related_terms if rel_term[0] == "is_a"] #Identify parent HPO term(s)
	if next_terms==["HP:0000118"]: #Level 0 term to indicate that the top-level term has been reached
		term_name=term.name
		top_level_list.append(term_name) #Add top-level term to list
	else:
		for next_term in next_terms: #Iterate through parent terms
			if next_term!="HP:0000118": #Level 0 term to indicate that the top-level term has been reached
				top_level_list=getNextTerms(next_term,top_level_list) #Recursion--move up 1 level in ontology
	return top_level_list #Return list to main method


#Import libraries
import csv
import sys
from orangecontrib.bio import ontology

#Open input phenotype file based on name provided at command line
batch=sys.argv[1]
infile=open(batch+".txt",'rU')
in_lines=infile.readlines()[1:]

#Initialize output file with header
outwriter = csv.writer(file(batch+"_hpo.txt",'w'),delimiter='\t')
header=["CNV","sample","input_phenotypes","input_hpo_terms","top_level_phenotypes"]
outwriter.writerow(header)

#Import human phenotype ontology file into an ontology data structure
ontology=ontology.OBOOntology("hp.obo")

#Iterate through all samples in input file
for line in in_lines:
	line=line.strip().split("\t")
	cnv=line[0] #Pathogenic CNV name
	sample=line[1] #Decipher sample ID
	phenotype_list=line[9].strip().split(", ") #Comma-separated list of phenotypes
	
	#Initialize lists of phenotypes and HPO terms for output
	phenotype_list_new=[]	
	hpo_term_list=[]
	top_term_list=[]

	for phenotype in phenotype_list: #Iterate through each phenotype for a sample
		phenotype=phenotype.strip('"')
		phenotype_list_new.append(phenotype)

		#Skip phenotype modifiers not prsent for HPO
		if phenotype=='' or phenotype=="mild" or phenotype=="moderate" or phenotype=="severe" or phenotype=="borderline" or phenotype=="profound" or phenotype=="Familial predisposition" or phenotype=="progressive":
			continue

		else:
			try: #Look up term in HPO data structure
				hpo_term=ontology.term_by_name(phenotype)
				hpo_term_list.append(hpo_term.id)
				
				top_level=[]
				top_level=getNextTerms(hpo_term,top_level) #Call method to get top-level HPO term for each phenotype
				for top_term in top_level:
					top_term_list.append(top_term)

			except ValueError: #Print error if term not in HPO
				print "Unidentified phenotype in sample "+sample+": "+phenotype
				continue

	#Output results to file
	top_term_list=list(set(top_term_list)) #Remove duplicate top terms for final output
	phenotype_string=','.join(phenotype_list_new)
	hpo_string=','.join(hpo_term_list)
	outrow=[cnv,sample,phenotype_string,hpo_string]+top_term_list
	outwriter.writerow(outrow)
