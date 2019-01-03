import re
import os
import argparse
#import csv


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


def extractAnnotation(annotationFile, output, SNP = None, region = None, extract = None, configureFile = None, bgzipOut = False):

	# save header as txt
	cmd_header = "tabix {0} -H".format(annotationFile) + ">" + "header.txt"
        os.system(cmd_header)

	with open("header.txt", 'r') as f:
                for line in f:
                        if "#" in line:
                                header = line.split("\t")
	

	#if given SNP chr:pos   
        if SNP != None:
                snp_split = SNP.split(":")
                cmd_snp = "tabix {0} {1}-{2} -h".format(annotationFile, SNP, snp_split[1])+ ">" + "tabix_result.txt"
		os.system(cmd_snp)
		

	# if given a region chr:posStart-posEnd
	if region != None:
		cmd_region = "tabix {0} {1} -h".format(annotationFile, region)+ ">" + "tabix_result.txt"
        	os.system(cmd_region)

	#if given bed file
	if extract != None:
		cmd_extract = "tabix {0} -B {1} -h".format(annotationFile, extract)+ ">" + "tabix_result.txt"
                os.system(cmd_extract)
	

	# default columns	
	column = ["# [1]CHROM", "[2]POS", "[3]REF", "[4]ALT", "[5]FILTER","ANNOVAR_refseq_Effect", "ANNOVAR_refseq_Transcript_ID", "ANNOVAR_refseq_Gene_ID", "ANNOVAR_refseq_Closest_gene(intergenic_only)", "ANNOVAR_refseq_HGVSc", "ANNOVAR_refseq_HGVSp", "ANNOVAR_refseq_Exon_Rank", "ANNOVAR_refseq_summary", "SnpEff_refseq_Effect", "SnpEff_refseq_Effect_impact", "SnpEff_refseq_Sequence_feature", "SnpEff_refseq_Sequence_feature_impact", "SnpEff_refseq_Transcript_ID", "SnpEff_refseq_Transcript_biotype", "SnpEff_refseq_Gene_name", "SnpEff_refseq_Gene_ID", "SnpEff_refseq_HGVSc", "SnpEff_refseq_HGVSp", "SnpEff_refseq_Protein_position/Protein_len", "SnpEff_refseq_CDS_position/CDS_len", "SnpEff_refseq_cDNA_position/cDNA_len", "SnpEff_refseq_Exon_or_intron_rank/total", "SnpEff_refseq_Distance_to_feature", "SnpEff_refseq_Warnings", "SnpEff_refseq_LOF/NMD", "SnpEff_refseq_LOF/NMD_gene_name", "SnpEff_refseq_LOF/NMD_gene_ID", "SnpEff_refseq_LOF/NMD_num_transcripts_affected", "SnpEff_refseq_LOF/NMD_percent_transcripts_affected", "SnpEff_refseq_TF_binding_effect", "SnpEff_refseq_TF_name", "SnpEff_refseq_TF_ID", "SnpEff_refseq_summary", "clinvar_rs", "clinvar_clnsig", "clinvar_trait", "clinvar_golden_stars", "GWAS_catalog_rs", "GWAS_catalog_trait", "GWAS_catalog_pubmedid", "1000Gp3_EUR_AF", "1000Gp3_AFR_AF", "ExAC_AFR_AF", "ExAC_NFE_AF", "CADD_phred", "SIFT_pred", "Eigen-phred", "Polyphen2_HVAR_pred"]
	
	# if specified a configure file
	if configureFile != None:
		find_col = open(configureFile, "r").read().rstrip("\n")
		column = find_col.split("\n") # Can change delimiter here for read-in txt

        column_int = []
        for item in column:
                if item in header:
                        column_int.append(header.index(item) + 1 )
	
	find_col_num = ",".join(map(str,column_int))

#	saved as csv file
#	temp_tabix = csv.writer(open("/udd/reliu/WGSA_extraction/temp.csv", "w"))
#	temp_tabix.writerows(csv.reader(open("/udd/reliu/WGSA_extraction/tabix_result.txt", "rb"),delimiter = '\t'))
	
	# extract columns and save as a txt file
	cmd_col = "cut -f {0} tabix_result.txt".format(find_col_num) + ">" + output
	os.system(cmd_col)

	if bgzipOut == True:
		cmd_bgzip = "bgzip -c {0} > {1}.bgz".format(output, output)
		os.system(cmd_bgzip)
		os.remove(output)
	
	# Convert txt to csv
#	txt_file = r"/udd/reliu/WGSA_extraction/txt_output.txt"
#	csv_file = r"/udd/reliu/WGSA_extraction/csv_output.csv"

#	in_txt = csv.reader(open(txt_file, "rb"), delimiter = '\t')
#	out_csv = csv.writer(open(csv_file, "w"))
		
#	out_csv.writerows(in_txt) 
		
		
	os.remove("header.txt")
	os.remove("tabix_result.txt")
			
	
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description = "WGSA extraction")
	subparser = parser.add_subparsers(help = "help")
	run_pipeline = subparser.add_parser("run", help = "run WGSA extraction")
	run_pipeline.set_defaults(which = "run")
	run_pipeline.add_argument("--annotationFile", required = True, help = "Annotation file")
	run_pipeline.add_argument("--region",required = False, help = "Region to bed extracted- chr:posStart-posEnd")
	run_pipeline.add_argument("--SNP", required = False, help = "SNP to be extracted- chr:pos")
	run_pipeline.add_argument("--extract", required = False, help = "A bed file for variants to be extracted")
	run_pipeline.add_argument("--configureFile", required = False, help = "Configure file")
	run_pipeline.add_argument("--out", required = True, help = "Output a txt file")
	run_pipeline.add_argument("--bgzipOut", required = False, action = "store_true", help = "Output a bgzipped txt file")
	
	args = vars(parser.parse_args())
	if args["which"] == "run":
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

	extractAnnotation(annoFile, output, SNP, region, extract, configureFile, bgzipOut)

	
