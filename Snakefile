rule setup:
	output:
		"/home/julius/Documents/results/cr/1000g/1000G_EUR_chr1.vcf.gz",
		"/home/julius/Documents/results/cr/1000g/genetic_map_GRCh37_fixed_chr1.txt",
		"/home/julius/Documents/results/cr/Homo_sapiens.GRCh37.genes"
	shell:
		"""
		# Install plink yourself if don't have it yet.
		# Install bcftools yourself if don't have it yet.

		cd ~/Documents/results/cr/
		# download ldsc:
		git clone https://github.com/bulik/ldsc.git
		cd ldsc; conda env create --file environment.yml
		# some packages need manual fixing for some reason:
		conda activate ldsc
		pip install cython
		pip install -Iv pandas==0.21

		# extract genes from the full gff
		wget https://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
		gunzip Homo_sapiens.GRCh37.87.gff3.gz -c | awk '$0~/^##/{next} $3=="gene"{print}' > Homo_sapiens.GRCh37.genes
		
		cd ~/Documents/results/cr/1000g/
		# Download 1000G EUR files from https://zenodo.org/record/6614170
		# (503 samples x 1.8 M snps)
		# wget https://zenodo.org/record/6614170/files/1000G_EUR.bed?download=1 -O 1000G_EUR.bed
		# wget https://zenodo.org/record/6614170/files/1000G_EUR.bim?download=1 -O 1000G_EUR.bim
		wget https://zenodo.org/record/6614170/files/1000G_EUR.fam?download=1 -O 1000G_EUR.fam

		# Download 1000G chr1 vcf:
		wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O 1000G_chr1.vcf.gz

		# also the ldscores:
		wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
		tar -jxvf eur_w_ld_chr.tar.bz2
		rm eur_w_ld_chr.tar.bz2
		gunzip genetic_map_GRCh37_chr1.txt.gz
		cut -f 2,3,4 genetic_map_GRCh37_chr1.txt > genetic_map_GRCh37_fixed_chr1.txt

		# extract Chr1 EURs and SNPs in ldscore files
		# (though SNP filter is not necessary as we recalculate our LD scores anyway)
		cut -f2 1000G_EUR.fam > EUR.txt
		bcftools view -S EUR.txt -i ID==@eur_w_ld_chr/w_hm3.snplist 1000G_chr1.vcf.gz | bcftools filter -i 'EUR_AF>0.05 && EUR_AF<0.95' | gzip > 1000G_EUR_chr1.vcf.gz

		"""

rule simulategenos:
	# Using sim1000G package.
	input:
		"/home/julius/Documents/results/cr/1000g/1000G_EUR_chr1.vcf.gz",
		"/home/julius/Documents/results/cr/1000g/genetic_map_GRCh37_fixed_chr1.txt"
	output:
		"/home/julius/Documents/results/cr/1000g/simulated_{n}.bed",
		"/home/julius/Documents/results/cr/1000g/simulated_{n}R.raw"
	params:
		bedstem="1000g/simulated",
		n=lambda wc: wc.get("n")
	shell:
		"""
		# Simulate ~90k snps from chr1 following 1kg EUR structure, in blocks,
		# based on the provided VCF:
		Rscript simulate-geno.R {params.n}
		# will create 1000g/tmp/simulated*map,ped.
		# Combine them into one bedfile
		# and attach cM values from 1kg EUR:
		cd ~/Documents/results/cr/
		paste <(ls 1000g/tmp/*ped) <(ls 1000g/tmp/*map) > 1000g/tmp/filelist.txt
		~/soft/plink1.9/plink --merge-list 1000g/tmp/filelist.txt \
			--cm-map 1000g/genetic_map_GRCh37_fixed_chr1.txt 1 \
			--make-bed --out {params.bedstem}_{params.n}
		# prepare a copy for loading to R
		~/soft/plink1.9/plink --bfile {params.bedstem}_{params.n} \
			--recode A --out {params.bedstem}_{params.n}R
		rm 1000g/tmp/simulated*.ped
		rm 1000g/tmp/simulated*.map
		"""

rule simulatestudygenos:
	input:
		"/home/julius/Documents/results/cr/1000g/simulated_1000.bed",
		"/home/julius/Documents/results/cr/1000g/simulated_3000.bed",
		"/home/julius/Documents/results/cr/1000g/simulated_10000.bed"
		
rule simulatephenos:
	input:
		"/home/julius/Documents/results/cr/1000g/simulated.bed"
	output:
		"/home/julius/Documents/results/cr/simphenos/cont.pheno"
	script:
		"./simulate-pheno.R"


# output:	~/Documents/results/cr/1000g/simulated_ldscores.ld.ldscore.gz
#		# Calculate LD scores
#		conda activate ldsc
#		./ldsc/ldsc.py --bfile {params.bedstem} \
#			--l2 \
#			--ld-wind-cm 1\
#			--out 1000g/simulated_ldscores
		
