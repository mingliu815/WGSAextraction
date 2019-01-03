import re
import os
import argparse
import subprocess
import pandas as pd
from pandas import DataFrame
from os.path import basename
from os.path import splitext


###### A python script for extraction of WGSA annotation files

##### Input files
## @ Indexed and sorted bgzipped WGSA annnotation files

##### Command options
## @ annotationFiles file1.bgz, file2.bgz
## @ SNP chr:pos (SNP to be extracted, optional)
## @ region chr:posStart-posEnd (region to bed extracted, optional)
## @ extract positionsNeeded.bed (or can specify a bed file for variants to be extracted, optional)
## @ configureFile columnsNeeded.txt (Optional)
## @ out output.csv (output.csv.gz)
## @ bgzipOut (if this is specified, the output file is bgzipped and indexed)


def extractAnnotation(input_num, inputfiles, output, SNP = None, region = None, extract = None, configureFile = None, bgzipOut = False):
	inputfile_list = inputfiles.split(",")
	output_list = []
	for i in range(int(input_num)):
		annotationFile = inputfile_list[i]
		temp_out = splitext(basename(annotationFile))[0] 
		all_header = temp_out + "_header"
		tabix_result = temp_out + "_tabix_result"
		res_header = temp_out + "_res_header"
		res_body = temp_out + "_res_body"
		res_out = temp_out + "_tempOut"

	# save header as txt
		cmd_header = "tabix {0} CHROM:POS ".format(annotationFile) + ">" + all_header
        	os.system(cmd_header)

		with open(all_header, 'r') as f:
                	for line in f:
                		 header = line.split("\t")

	#if given SNP chr:pos   
        	if SNP != None:
                	snp_split = SNP.split(":")
                	cmd_snp = "/udd/reliu/bin/htslib/tabix {0} {1}-{2} -h".format(annotationFile, SNP, snp_split[1])+ ">" + tabix_result
			os.system(cmd_snp)
		

	# if given a region chr:posStart-posEnd
		if region != None:
			cmd_region = "/udd/reliu/bin/htslib/tabix {0} {1} -h".format(annotationFile, region)+ ">" + tabix_result
        		os.system(cmd_region)

	#if given bed file
		if extract != None:
			bedtolist = temp_out + "_CmdBedtoList.txt"
			cmd_toList = "Rscript /udd/reliu/WGSA_extraction/BedToList.R {0} ".format(extract) + ">" + bedtolist
			os.system(cmd_toList)
			cmd_extract = "xargs -a {0} -I {{}} /udd/reliu/bin/htslib/tabix {1} {{}} > {2}".format(bedtolist, annotationFile, tabix_result)
                	os.system(cmd_extract)
			os.remove(bedtolist)

	# default columns	
		column = ["# [1]CHROM", "[2]POS", "[3]REF", "[4]ALT", "[5]FILTER","ANNOVAR_refseq_Effect", "ANNOVAR_refseq_Transcript_ID", "ANNOVAR_refseq_Gene_ID", "ANNOVAR_refseq_Closest_gene(intergenic_only)", "ANNOVAR_refseq_HGVSc", "ANNOVAR_refseq_HGVSp", "ANNOVAR_refseq_Exon_Rank", "ANNOVAR_refseq_summary", "SnpEff_refseq_Effect", "SnpEff_refseq_Effect_impact", "SnpEff_refseq_Sequence_feature", "SnpEff_refseq_Sequence_feature_impact", "SnpEff_refseq_Transcript_ID", "SnpEff_refseq_Transcript_biotype", "SnpEff_refseq_Gene_name", "SnpEff_refseq_Gene_ID", "SnpEff_refseq_HGVSc", "SnpEff_refseq_HGVSp", "SnpEff_refseq_Protein_position/Protein_len", "SnpEff_refseq_CDS_position/CDS_len", "SnpEff_refseq_cDNA_position/cDNA_len", "SnpEff_refseq_Exon_or_intron_rank/total", "SnpEff_refseq_Distance_to_feature", "SnpEff_refseq_Warnings", "SnpEff_refseq_LOF/NMD", "SnpEff_refseq_LOF/NMD_gene_name", "SnpEff_refseq_LOF/NMD_gene_ID", "SnpEff_refseq_LOF/NMD_num_transcripts_affected", "SnpEff_refseq_LOF/NMD_percent_transcripts_affected", "SnpEff_refseq_TF_binding_effect", "SnpEff_refseq_TF_name", "SnpEff_refseq_TF_ID", "SnpEff_refseq_summary", "clinvar_rs", "clinvar_clnsig", "clinvar_trait", "clinvar_golden_stars", "GWAS_catalog_rs", "GWAS_catalog_trait", "GWAS_catalog_pubmedid", "1000Gp3_EUR_AF", "1000Gp3_AFR_AF", "ExAC_AFR_AF", "ExAC_NFE_AF", "CADD_phred", "SIFT_pred", "Eigen-phred", "Polyphen2_HVAR_pred"]
	
	# if specified a configure file
		if configureFile != None:
			find_col = open(configureFile, "r").read().rstrip("\n")
			column = find_col.split("\n") # Can change delimiter here for read-in txt

        	column_int = []
		column_specific = []
        	for item in column:
                	if item in header:
				column_specific.append(item)
                        	column_int.append(header.index(item) + 1 )
	
		find_col_num = ",".join(map(str,column_int))

#	saved as csv file
#	temp_tabix = csv.writer(open("/udd/reliu/WGSA_extraction/temp.csv", "w"))
#	temp_tabix.writerows(csv.reader(open("/udd/reliu/WGSA_extraction/tabix_result.txt", "rb"),delimiter = '\t'))
	
	# extract columns and save as a txt file
		cmd_col = "cut -f {0} {1}".format(find_col_num, tabix_result) + ">" + res_body
		os.system(cmd_col)
	
		header_out = "\t".join(column_specific)
		f = open(res_header,'w')
		f.write(header_out + "\n")
		f.close()
		cmd_cat = "cat {0} {1} > {2}".format(res_header, res_body, res_out)
		os.system(cmd_cat)

#		if bgzipOut == True:
#			cmd_bgzip = "bgzip -c {0} > {1}.bgz".format(output, output)
#			os.system(cmd_bgzip)
#			os.remove(output)
		
		
		os.remove(all_header)
		os.remove(tabix_result)
		os.remove(res_header)
		os.remove(res_body)			
		output_list.append(res_out)

	if len(output_list) == 1:
		cmd_out = "cp {0} {1}".format(output_list[0], output)
		os.system(cmd_out)

	elif len(output_list) == 2:
		df1 = pd.read_csv(output_list[0], sep="\t")
		df2 = pd.read_csv(output_list[1], sep="\t")

		common_col = list(set(list(df1))&set((list(df2))))

		df1["sorted_id"] = df1["CHROM"].astype(str)+":"+df1["POS"].astype(str)+":"+df1["REF"].astype(str)+":"+df1["ALT"].astype(str)
		df2["sorted_id"] = df2["CHROM"].astype(str)+":"+df2["POS"].astype(str)+":"+df2["REF"].astype(str)+":"+df2["ALT"].astype(str)

		df1.drop(common_col, inplace=True, axis=1)

#dfout = df1.set_index("sorted_id").join(df2.set_index("sorted_id"))
		dfout = df1.merge(df2, on='sorted_id')

		dfout.drop("sorted_id", inplace=True, axis=1)

#dfout.to_csv('/udd/reliu/WGSA_extraction/newOutput.txt', sep='\t')
		dfout.to_csv(output, sep='\t', index=False)

	if bgzipOut == True:
		cmd_bgzip = "bgzip -c {0} > {1}.bgz".format(output, output)
                os.system(cmd_bgzip)
                os.remove(output)

	for item in output_list:
		os.remove(item)	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description = "WGSA extraction")
#	subparser = parser.add_subparsers(help = "help")
#	run_pipeline = subparser.add_parser("run", help = "run WGSA extraction")
#	run_pipeline.set_defaults(which = "run")
	parser.add_argument("--fileNum", required = True, help = "Number of annotation files")
	parser.add_argument("--annotationFile", required = True, help = "Path of annotation files, split by ','")
	parser.add_argument("--region",required = False, help = "Region to bed extracted- chr:posStart-posEnd")
	parser.add_argument("--SNP", required = False, help = "SNP to be extracted- chr:pos")
	parser.add_argument("--extract", required = False, help = "A bed file for variants to be extracted")
	parser.add_argument("--configureFile", required = False, help = "Configure file")
	parser.add_argument("--out", required = True, help = "Output a txt file")
	parser.add_argument("--bgzipOut", required = False, action = "store_true", help = "Output a bgzipped txt file")
	
	args = vars(parser.parse_args())
#	if args["which"] == "run":
	inputNum = args["fileNum"]
	annoFile = args["annotationFile"]
	region = args["region"]
	SNP = args["SNP"]
	extract = args["extract"]
	configureFile = args["configureFile"]
	output = args["out"]
	bgzipOut = False
	if region != None:
		region = region
	if SNP != None:
		SNP = SNP
	if extract != None:
		extract = extract
	if configureFile != None:
                configureFile = configureFile
	if parser.parse_args().bgzipOut:
		bgzipOut = True

	extractAnnotation(inputNum, annoFile, output, SNP, region, extract, configureFile, bgzipOut)
