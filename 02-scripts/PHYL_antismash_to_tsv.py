'''
Authors: Vittorio Tracanna, Jasmin Frangenberg
Language: Python
Objective:
	Extract the following information from multiple antiSMASH results:
	- Sample name (which is present in the name of the summary GBK file or the name of the antiSMASH folder)
	- Contig ID of BGC region (which is present in the filename of the individual BGC region gbk files)
	- BGC type
	- Is BGC region on contig edge yes/no?
	- Length of BGC core region and whole BGC region
	- knownclusterblast result (maybe even combined with your refinement script)
	- How many CDSs are in the BGC region
Output: one TSV file
'''
print("Script for extracting results from antiSMASH output")

import os
import re
import sys

# Check if BioPython is installed. Otherwise offer the user to automatically install it
try:
	from Bio import SeqIO # Try to import the BioPython module
except ImportError: # If BioPython is not installed, print these options to the user
	print(
	"\nThis script requires the module BioPython. Looks like it is not installed. I can",
	"install it for you in your current Conda environment (~17MB of disk space).",
	"If you want to create a separate Conda environment first, please do that and",
	"then run this script again. For example:",
	"\t> conda create -n biopythonEnv python biopython",
	"\t> conda activate biopythonEnv",
	"\t> python " + sys.argv[0],
	"If you don't use Conda you can get BioPython via 'pip install biopython'.\n",
	sep = "\n")

	ans = ""
	while ans != "y": # Wait for user decision
		ans = input("Do you want me to install BioPython? (y/n) ").strip()
		if ans == "y": # If yes: Start BioPython installation
			import subprocess
			try:
				subprocess.check_call(["conda", "install", "biopython"])
				try:
					from Bio import SeqIO
				except ModuleNotFoundError: # Tell the user if there was an error during installation and exit.
					print("BioPython installation aborted.")
					exit()
			except subprocess.CalledProcessError: # BioPython installation will abort if the current Python version is not compatible to BioPython. Let the user know and exit.
				print(
				"Arrrg. Your Python version does not fit (it is probably too new). You can make",
				"sure to have compatible versions of Python and BioPython by installing them",
				"in a new Conda environment using the example commands from above.",
				"Program aborted.",
				sep = "\n")
				exit()
		elif ans == "n": # If user does not want us to install BioPython: Exit
			print("Program aborted.")
			exit()

# Function extract_gbk_info:
# - Iterates over folders with antiSMASH output
# - Opens each GBK and parses the information
# - Extracts the knownclusterblast output from the antiSMASH folder
# - Stores everything into a dictionary that can then be written on disk and easily opened in Excel, Calc or R
def extract_gbk_info(as_path, dirs, contigs):
	as_results = {}
	for antismash_output_folder in dirs:
		error_ids = []
		as_dir_path = '{}{}/'.format(as_path, antismash_output_folder) # Path without GBK file
		print("  - " + antismash_output_folder)
		gbks = [x for x in os.listdir(as_dir_path) if x.endswith('.gbk') and 'region' in x]

		for gbk in gbks: # Go through all the antiSMASH output gbk files
			as_file_path = '{}{}'.format(as_dir_path, gbk) # Full path to current gbk

			# Check if the antiSMASH output looks as expected
			record_num = 0 # Enumerate GenBank records. There should only be one per file (= one BGC per antiSMASH output GBK file)
			for record in SeqIO.parse(open(as_file_path), "genbank"): # Parse antiSMASH output file
				record_num += 1
				if record_num > 1:
					print("Here is something fishy. Is there more than 1 BGC in this GBK file?",
					as_file_path,
					sep="\n")
					print("You can continue or abort. Continuing might lead to unexpected results as this script was not tested for multiple GenBank records in antiSMASH output.\n")
					if input("Continue? (y/n) ").strip() == "n":
						print("Program aborted.")
						exit()

				# Start parsing and assigning variables
				Sample_name = antismash_output_folder
				Contig_id = record.annotations["accessions"][0]

				# Go through the features to look at cand_cluster and CDS
				cand_cluster_count = 0
				CDS_count = 0 # Count the open reading frames of the cluster
				for feature in record.features:

					# Use the first cand_cluster feature from the GBK file to extract all the infos; the first one holds a summary of all properties found by antiSMASH:
					# - length (BGC_length)
					# - number of all products (num_BGC_core_genes)
					# - all products (BGC_core_genes). Store them in a string (comma-separatedly)
					# Write a single line for the BGC region in the TSV file
					if feature.type == "cand_cluster" and cand_cluster_count == 0: # we dont need cand_cluster feature anymore. replace with protocluster
						cand_cluster_count += 1
						BGC_length = str(len(feature.extract(record.seq)))
						num_BGC_core_genes = str(len(feature.qualifiers["product"]))
						BGC_core_genes = feature.qualifiers["product"]
						for i in range(len(BGC_core_genes)):
							BGC_core_genes[i] = BGC_core_genes[i][0].upper() + BGC_core_genes[i][1:] # Make first letters uppercase, e.g. lassopeptide -> Lassopeptide

					# Count functional CDSs (no pseudogenes)
					elif feature.type == "CDS" and "translation" in feature.qualifiers.keys(): # Make sure not to count pseudogenes (which would have no "translation tag")
						CDS_count += 1

			infile = '{}knownclusterblast/{}_c1.txt'.format(as_dir_path, record.id)
			clust_IDs, clust_annotations, blast_num, blast_identity_averages, blast_score_cum, blast_score_averages, blast_coverage_averages, blast_cds_annotations = parse_cluster(infile)
			BGC_pos = record.annotations["structured_comment"]
			contig_edge, error_id = edge_position(contigs, antismash_output_folder, record.id, BGC_pos) # Get the side of the contig edge (left/right/both/none) #### Buth this into as_results and TSV output
			error_ids.append(error_id) # Contig IDs for which no prokka output was found

			# Store the all the values for current GBK in a list
			as_results[gbk] = [Sample_name,
								Contig_id,
								num_BGC_core_genes,
								BGC_core_genes,
								contig_edge,
								BGC_length,
								str(CDS_count),
								clust_IDs, # IDs and annotations of BGC BLAST matches
								clust_annotations,
								blast_num,
								blast_identity_averages,
								blast_score_cum,
								blast_score_averages,
								blast_coverage_averages,
								blast_cds_annotations]
		
		# Print contig IDs which are not found in the prokka output to standard output. Print 3 contig IDs per line.
		error_ids = [error_id for error_id in error_ids if error_id] # Remove empty strings
		if len(error_ids) == 1:
			print("      Warning: The following contig was not found in prokka fna files.")
			print("               No contig edge inference possible for:")
			print("               " + error_ids[0])
		elif len(error_ids) > 1:
			print("      Warning: The following contigs were not found in prokka fna files.")
			print("               No contig edge inference possible for:")
			while len(error_ids) > 0:
				try:
					print("               " + error_ids[0] + ", " + error_ids[1] + ", " + error_ids[2], end="")
					if len(error_ids) > 3:
						print(",") # Put comma at the end of line
						error_ids = error_ids[3:]
					else:
						print()
						error_ids = []
				except IndexError: # If only two IDs are left
					try:
						print("               " + error_ids[0] + ", " + error_ids[1])
						error_ids = []
					except IndexError: # If only one ID is left
						print("               " + error_ids[0])
						error_ids = []

	return as_results

# Function: Get contig lengths from prokka contigs FASTA files
def get_contig_lengths(pr_path, dirs):
	d = {} # Dictionary of contig IDs and their lengths, grouped in a dictionary with sample IDs as keys
	for dir in dirs:
		d[dir] = {}
		pr_dir_path = '{}{}/'.format(pr_path, dir) # Path without fna file
		fnas = [x for x in os.listdir(pr_dir_path) if x.endswith('.fna')]

		# Parse fna files
		for fna in fnas:
			pr_file_path = '{}{}'.format(pr_dir_path, fna) # Full path to current fasta file
			fasta = SeqIO.parse(open(pr_file_path), "fasta")
			for contig in fasta:
				d[dir][contig.id] = len(contig.seq)
	return d

# Function: Determine on which side of the contig the BGC stops. Region of 50 base pairs on each end counts as truncated
def edge_position(contigs, sample, id, bgc_pos):
	bgc_start = int(bgc_pos["antiSMASH-Data"]["Orig. start"]) # BGC start position in contig
	bgc_end = int(bgc_pos["antiSMASH-Data"]["Orig. end"]) # BGC stop position in contig
	
	try:
		if bgc_start <= 50 and bgc_end >= contigs[sample][id] - 50: # Contig edge on both sides
			side = "left, right"
		elif bgc_start > 50 and bgc_end < contigs[sample][id] - 50: # No contig edge
			side = "no"
		elif bgc_start <= 50 and bgc_end < contigs[sample][id] - 50: # Contig edge on the left
			side = "left"
		elif bgc_start > 50 and bgc_end >= contigs[sample][id] - 50: # Contig edge on the right
			side = "right"
		error = ""
	except KeyError: # Return contig IDs for which no entry in the prokka file was found
		side = "NA"
		error = id

	return side, error

# Function: Read in the antiSMASH txt file and split it into the records (divided by ">>")
def file_to_list(file_path):
	hits = open(file_path).read()

	# Divide the hits into a list of lists:
	# First list level: All the records divided by ">>"
		# Second list level: Metadata at the top of each record, genes with annotations, significant BLAST hits
		# Final list looks like this: [[Metadata1, genes with annotation1, BLAST hits1], [Metadata2, genes with annotation2, BLAST hits2], [...]]
	hits = hits.split("\n\n>>\n")[1:] # Split by record
	for i in range(len(hits)): # Split each record into metadata, annotations and significant BLAST hits
		hits[i] = hits[i].split("\nTable of genes, locations, strands and annotations of subject cluster:\n")
		metadata = hits[i][0].rstrip().split("\n")

		annotation_data = hits[i][1].split("\nTable of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):\n")
		annotations = annotation_data[0].rstrip().split("\n")
		annotations = [gene.split("\t") for gene in annotations]

		blast_hits = annotation_data[1].rstrip().split("\n")
		blast_hits = [hit.split("\t") for hit in blast_hits]
		hits[i] = [metadata, annotations, blast_hits]
	return[hits]

# Function: Parse the antiSMASH KnownclusterBLAST txt files
def parse_cluster(file_path): # Function for calculating the average BLAST hits identity value and extracting datails from the txt files

	# Get a list of all the records in the antiSMASH txt file
	records = file_to_list(file_path)

	# Extract and calculate the results from that list
	## (Known)clusterblast ID
	clust_IDs = [re.search("\d+\. (.*)", metaline[0][0]).group(1) for metaline in records[0]]

	## (Known)clusterblast annotations
	clust_annotations = [re.search("Source: (.*)", metaline[0][1]).group(1) for metaline in records[0]]
	for j in range(len(clust_annotations)):  # Make annotations uppercase
		clust_annotations[j] = clust_annotations[j][0].upper() + clust_annotations[j][1:]

	## Cumulative BLAST score
	blast_score_cum = [re.search("Cumulative BLAST score: (\d+)", metaline[0][4]).group(1) for metaline in records[0]]

	## Average BLAST identity to (Known)clusterblast genes
	blast_num = []
	blast_identity_averages = []
	blast_score_averages = []
	blast_coverage_averages = []
	blast_cds_annotations = []

	for record in records[0]:

		# Store (Known)clusterblast IDs and corresponding annotations in dictionary (if present)
		try:
			gene_match_ids = {line[0]: line[5] for line in record[1]} # Dictionary of (Known)clusterblast gene IDs and annotations
		except IndexError:
			gene_match_ids = {}

		# Extract annotation IDs, identity values, BLAST scores and subject coverage from BLAST hits table
		blast_hit_ids = [line[1] for line in record[2]]
		blast_hit_identities = [float(line[2]) for line in record[2]]
		blast_hit_scores = [float(line[3]) for line in record[2]]
		blast_hit_coverages = [float(line[4]) for line in record[2]]

		# Get number of BLAST hits per record
		blast_num.append(str(len(record[2])))

		# Get (Known)clusterblast cds_annotations
		blast_cds_annotations.append([gene_match_ids[hit_id] for hit_id in blast_hit_ids if hit_id in gene_match_ids.keys()])

		# Calculate avarage identity score
		blast_identity_averages.append(str(round(sum(blast_hit_identities)/len(blast_hit_identities), 2)))
		blast_score_averages.append(str(round(sum(blast_hit_scores)/len(blast_hit_scores), 2)))
		blast_coverage_averages.append(str(round(sum(blast_hit_coverages)/len(blast_hit_coverages), 2)))

	return clust_IDs, clust_annotations, blast_num, blast_identity_averages, blast_score_cum, blast_score_averages, blast_coverage_averages, blast_cds_annotations

# Check if input/output folders are given on the command line
if len(sys.argv) != 4:
	print(
		"\nYou need to provide the prokka and antiSMASH directories as well as an output file name.", #### take care of character limit (80)
		"You can enter them here or cancel (Ctrl + c) and run the script again like:\n",
		"\tpython " + sys.argv[0] + " [antiSMASH directory] [output file]\n",
		sep="\n")
	try:
		as_path = ""
		while as_path == "":
			as_path = input("antiSMASH directory: ").strip()
		pr_path = ""
		while pr_path == "":
			pr_path = input("Prokka directory: ").strip()
		outpath = input("Output file: ").strip() ##### Take care of "" as output (also by loop?)
	except KeyboardInterrupt:
		exit("\nProgram aborted. Please execute the script again.")
else:
	pr_path = sys.argv[1] # Path to the folder where all the prokka contigs are
	as_path = sys.argv[2] # Path to the folder where all the antiSMASH results are
	outpath = sys.argv[3] # Path to the TSV file

# Store all the prokka contig lengths in dictionary contigs
if not pr_path.endswith("/"): # Make sure the input path ends with a slash (for later directory parsing)
	pr_path += "/"
try:
	pr_dirs = [d for d in os.listdir(pr_path) if os.path.isdir('{}{}'.format(pr_path, d))] # List of all the directories in the base path
	print("\nFound these prokka directories:")
	for dir in pr_dirs:
		print("  " + dir)
except FileNotFoundError:
	exit("No such directory '" + pr_path + "' found. Please check your directory name and try again.")

# Store all the antiSMASH results in dictionary as_results
if not as_path.endswith("/"): # Make sure the input path ends with a slash (for later directory parsing)
	as_path += "/"
try:
	as_dirs = [d for d in os.listdir(as_path) if os.path.isdir('{}{}'.format(as_path, d))] # List of all the directories in the base path
	print("\nFound these antiSMASH directories:")
	for dir in as_dirs:
		print("  " + dir)
except FileNotFoundError:
	exit("No such directory '" + as_path + "' found. Please check your directory name and try again.")

# Make sure all the prokka and antiSMASH directories match each other. If not: Don't parse the missing directories.
dirs = []
for pr_dir in pr_dirs:
	if pr_dir in as_dirs:
		dirs.append(pr_dir)

# Get dictionary of prokka contigs and their respective lengths
contigs = get_contig_lengths(pr_path, dirs)

# Dictionary with all the results to be written to disk
print("\nParsing these samples:")
as_results = extract_gbk_info(as_path, dirs, contigs)

# Write the antiSMASH results into TSV (tab-separated values) file
with open(outpath, "w", newline='') as outfile:
	fieldnames = ['Sample_name', 'Contig_id', 'Num_BGC_core_genes', 'BGC_core_genes', 'BGC_at_contig_edge', 'BGC_length', 'CDS_count', 'Knownclusterblast_ID', 'Knownclusterblast_annotation', 'BLAST_hits', 'Knownclusterblast_identity_average', "BLAST_score_cumulative", "BLAST_score_average", "BLAST_coverage_average", 'BLAST_annotations']
	outfile.write("\t".join(fieldnames) + "\n") # Write column names

	for sample in sorted(as_results.keys()): # Sort the samples not to have them in arbitrary order
		out_rows = len(as_results[sample][7]) # Determine if BLAST results were found - if yes: out_rows will be > 0
		if out_rows > 0: # If BLAST results are present: Write each result into a line in the TSV file
			for i in range(out_rows):
				for j in range(len(as_results[sample][3])): # If several products were found in the GBK file, write a line for each one
					outline1 = "\t".join(as_results[sample][:3])
					outline2 = as_results[sample][3][j]
					outline3 = "\t".join(as_results[sample][4:7])
					outline4 = as_results[sample][7][i]
					outline5 = as_results[sample][8][i]
					outline6 = as_results[sample][9][i]
					outline7 = as_results[sample][10][i]
					outline8 = as_results[sample][11][i]
					outline9 = as_results[sample][12][i]
					outline10 = as_results[sample][13][i]
					outline11 = ", ".join(as_results[sample][14][i])
					outfile.write("\t".join([outline1, outline2, outline3, outline4, outline5, outline6, outline7, outline8, outline9, outline10, outline11]) + "\n")
		else: # If no BLAST results are present: Write only the BGC into a line in the TSV file
			for j in range(len(as_results[sample][3])):	 # If several products were found in the GBK file, write a line for each one
				outline1 = "\t".join(as_results[sample][:3])
				outline2 = as_results[sample][3][j]
				outline3 = "\t".join(as_results[sample][4:7])
				outline4 = ""
				outline5 = ""
				outline6 = ""
				outline7 = ""
				outline8 = ""
				outline9 = ""
				outline10 = ""
				outline11 = ""
				outfile.write("\t".join([outline1, outline2, outline3, outline4, outline5, outline6, outline7, outline8, outline9, outline10, outline11]) + "\n")

print("\nFinished successfully. Your output is here: " + outpath)
