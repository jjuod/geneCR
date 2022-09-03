# 1. download 1000G EUR files from https://zenodo.org/record/6614170
# (503 samples x 1.8 M snps)
# also the ldscores:
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
# tar -jxvf eur_w_ld_chr.tar.bz2
# rm eur_w_ld_chr.tar.bz2

# 2. download ldsc:
# git clone https://github.com/bulik/ldsc.git
# cd ldsc; conda env create --file environment.yml
# (actually need to `conda activate ldsc` and then `pip install cython` and `pip install -Iv pandas==0.21` manually)

# 3. reduce the 1kg data a bit, to about 490k snps:
cd ~/Documents/results/cr/1000g
awk 'NR%2==1' eur_w_ld_chr/w_hm3.snplist > reduced_snplist.txt
~/soft/plink1.9/plink --bfile 1000G_EUR \
	--extract reduced_snplist.txt --maf 0.05 --geno 0 \
	--make-bed --out 1000G_sm

# 4. run the R script to simulate phenotype
# ./simulate-pheno.R
# analyze via plink
~/soft/plink1.9/plink --bfile 1000g/1000G_sm \
	--pheno simphenos/cont.pheno --assoc \
	--out simphenos/res_cont

# using t scores as z scores and fake alleles
awk 'NR==1{print "SNP", "A1", "A2", "N", "Z"; next} {print $2, "A1", "A2", $4, $8}' simphenos/res_cont.qassoc > simphenos/res_cont.ldready

# run LD score regr
./ldsc/ldsc.py --h2 simphenos/res_cont.ldready \
	--ref-ld-chr 1000g/eur_w_ld_chr/ \
	--w-ld-chr 1000g/eur_w_ld_chr/ \
	--out simphenos/ldsc

#... seems like 500 people isn't enough.

# ---------------
# OPTION B: sim1000G
# 1. still need the 1000G_EUR.fam file from https://zenodo.org/record/6614170
# then, download 1000G chr1 vcf from http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
# and extract EURs:
cd ~/Documents/results/cr/1000g
awk '{print $2}' 1000G_EUR.fam > EUR.txt
bcftools view -S EUR.txt -i ID==@eur_w_ld_chr/w_hm3.snplist 1000G_chr1.vcf.gz | bcftools filter -i 'EUR_AF>0.05 && EUR_AF<0.95' | gzip > 1000G_EUR_chr1.vcf.gz

# then use R to create 1000g/simulated.ped and .map
# then
~/soft/plink1.9/plink --file 1000g/simulated --make-bed simulated
rm 1000g/simulated.ped
rm 1000g/simulated.map

