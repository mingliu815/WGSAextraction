# WGSAextraction

WGSA (WGS annotator) is an annotation pipeline for human genome sequencing studies. More information is [here](https://sites.google.com/site/jpopgen/wgsa)

These are python scripts that extract information from a sorted and indexed bgzipped WGSA annotation files including both SNP and tissue-specific annotaion files.
Users can extract from a single file or two files (snp and tissued-specific) together by specifying the number of files. 
**If given two files at the same time, it will merge columns from two files into a final output.**  
WGSAextraction.py works for freeze2 and freeze3 data format and WGSAextraction\_freeze5.py works specifically for freeze5 data format.  

They accept information of CHR, start position and end position, and optionally a configuration file with selected columns. The output will be saved as a txt file (or a bgzipped txt file if specified) with all required information.

The configuration file must include **CHR, POS, REF and ALT** with other columns users specify.

The main extraction tool applied is tabix. 

### Command options

- "--annotationFile": give a path of an original sorted and indexed annotation file. **Required**
- "--SNP": extract a single SNP, with 'chr:pos' format. **Optional**
- "--region": extract a region of SNPs with 'chr:posStart-posEnd' format. **Optional**
- "--extract": give a bed file for variants to be extracted. Only first three fields are required: chrom, chromStart, and chromEnd. **Optional**
- "--configureFile": give a configuration file with selected columns in a txt file. If don't include this command, the output will include the default columns as below. **Optional**
- "--out": save output as a txt file. **Required**
- "--bgzipOut": save output as a bgzipped txt file (~.txt.bgz). **Optional** 

\*Default columns:
- CHROM
- POS
- REF
- ALT
- FILTER
- Annovar\_refseq\* (ANNOVAR\_refseq\_Effect, ANNOVAR\_refseq\_Transcript\_ID, ANNOVAR\_refseq\_Gene\_ID, ANNOVAR\_refseq\_Closest\_gene(intergenic\_only), ANNOVAR\_refseq\_HGVSc, ANNOVAR\_refseq\_HGVSp, ANNOVAR\_refseq\_Exon\_Rank, ANNOVAR\_refseq\_summary)
- SnpEff\_refseq\* (SnpEff\_refseq\_Effect, SnpEff\_refseq\_Effect\_impact, SnpEff\_refseq\_Sequence\_feature, SnpEff\_refseq\_Sequence\_feature\_impact, SnpEff\_refseq\_Transcript\_ID, SnpEff\_refseq\_Transcript\_biotype, SnpEff\_refseq\_Gene\_name, SnpEff\_refseq\_Gene\_ID, SnpEff\_refseq\_HGVSc, SnpEff\_refseq\_HGVSp, SnpEff\_refseq\_Protein\_position/Protein\_len, SnpEff\_refseq\_CDS\_position/CDS\_len, SnpEff\_refseq\_cDNA\_position/cDNA\_len, SnpEff\_refseq\_Exon\_or\_intron\_rank/total, SnpEff\_refseq\_Distance\_to\_feature, SnpEff\_refseq\_Warnings, SnpEff\_refseq\_LOF/NMD, SnpEff\_refseq\_LOF/NMD\_gene\_name, SnpEff\_refseq\_LOF/NMD\_gene\_ID, SnpEff\_refseq\_LOF/NMD\_num\_transcripts\_affected, SnpEff\_refseq\_LOF/NMD\_percent\_transcripts\_affected, SnpEff\_refseq\_TF\_binding\_effect, SnpEff\_refseq\_TF\_name, SnpEff\_refseq\_TF\_ID, SnpEff\_refseq\_summary)
- clinvar\* (clinvar\_rs, clinvar\_clnsig, clinvar\_trait, clinvar\_golden\_stars)
- GWAS\_catalog\* (GWAS\_catalog\_rs, GWAS\_catalog\_trait, GWAS\_catalog\_pubmedid)
- 1000Gp3\_EUR\_AF
- 1000Gp3\_AFR\_AF
- ExAC\_AFR\_AF
- ExAC\_NFE\_AF
- CADD\_phred
- SIFT\_pred
- Eigen-phred
- Polyphen2\_HVAR\_pred

\* Default columns are not suitable for freeze5 data format. Please make sure you have your own configuration file with matched column names when you use WGSAextraction\_freeze5.py 

### Examples

#### One annotation file

- If want to extract a single SNP with default columns:
	- Run `python WGSAextraction_freeze2.3.py run --annotationFile WGSAfile.bgz --SNP 1:99705 --out output.txt`
- If want to extract a region and specific columns needed in a bgzipped txt file:
	- Create a configuration file (sample as [columnsNeeded.txt])
	- Run `python WGSAextraction_freeze2.3.py run --annotationFile WGSAfile.bgz --region 1:1000-100000 --configureFile columnsNeeded.txt --out output.txt --bgzipOut` This will output a _output.txt.bgz_ file.
- If want to extract given a bed file and specific columns needed:
	- Create a bed file with only first three fields (sample as [positionsNeeded.bed])
	- Run `python /udd/reliu/WGSA_extraction/WGSAextraction_freeze2.3.py run --annotationFile WGSAfile.bgz --extract positionsNeeded.bed --configureFile columnsNeeded.txt --out output.txt`

#### SNP annotation and tissue-specifc annotation files together

- If want to do extraction from chromosome 1 for freeze5 given position needed and columns needed:
	- Configuration file for two files looks like [columnsNeeded\_snp.tissue.txt]).	
	- Run `python WGSAextraction_freeze5.py --fileNum 2 --annotationFile WGSAfile.snp.bgz,WGSAfile.tissue-specific.bgz --extract positionsNeeded.bed --configureFile columnsNeeded_snp.tissue.txt --out output.txt`
