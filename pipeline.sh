# 3-5. SIMULATE DATA AND APPLY LDSCORE
# OPTION A: directly using 1kg data:
if false
then
	# 3. reduce the 1kg data a bit, to about 490k snps:
	cd ~/Documents/results/cr/1000g
	awk 'NR%2==1' eur_w_ld_chr/w_hm3.snplist > reduced_snplist.txt
	~/soft/plink1.9/plink --bfile 1000G_EUR \
		--extract reduced_snplist.txt --maf 0.05 --geno 0 \
		--make-bed --out 1000G_sm
	
	# run the R script to simulate phenotype
	# ./simulate-pheno.R
	# 4. analyze via plink
	~/soft/plink1.9/plink --bfile 1000g/1000G_sm \
		--pheno simphenos/cont.pheno --assoc \
		--out simphenos/res_cont

	# using t scores as z scores and fake alleles
	awk 'NR==1{print "SNP", "A1", "A2", "N", "Z"; next} {print $2, "A1", "A2", $4, $8}' simphenos/res_cont.qassoc > simphenos/res_cont.ldready
	
	# 5. run LD score regr
	conda activate ldsc
	./ldsc/ldsc.py --h2 simphenos/res_cont.ldready \
		--ref-ld-chr 1000g/eur_w_ld_chr/ \
		--w-ld-chr 1000g/eur_w_ld_chr/ \
		--out simphenos/ldsc
fi
#... seems like 500 people isn't enough.

# OPTION B: sim1000G. See snakefile

# ---------------
# run the R script to simulate phenotype
# ./simulate-pheno.R

# analyze via plink
~/soft/plink1.9/plink --bfile 1000g/simulated \
	--pheno simphenos/cont.pheno --assoc \
	--out simphenos/res_cont
~/soft/plink1.9/plink --bfile 1000g/simulated \
	--pheno simphenos/bin.pheno --logistic beta \
	--out simphenos/res_bin


# using t scores as z scores and fake alleles
awk 'NR==1{print "SNP", "A1", "A2", "N", "Z"; next} {print $2, "A1", "A2", $4, $8}' simphenos/res_cont.qassoc > simphenos/res_cont.ldready
awk 'NR==1{print "SNP", "A1", "A2", "N", "Z"; next} {print $2, "A1", "A2", $6, $8}' simphenos/res_bin.assoc.logistic > simphenos/res_bin.ldready

# run LD score regr
conda activate ldsc
./ldsc/ldsc.py --h2 simphenos/res_cont.ldready \
	--ref-ld 1000g/simulated_ldscores \
	--w-ld 1000g/simulated_ldscores \
	--out simphenos/ldsc

./ldsc/ldsc.py --h2 simphenos/res_bin.ldready \
	--ref-ld 1000g/simulated_ldscores \
	--w-ld 1000g/simulated_ldscores \
	--samp-prev 0.10 --pop-prev 0.10 \
	--out simphenos/ldscb
# this kind of works and gives similar results.

# and this, blindly metaanalysing lin and log betas, gives worse results:
./ldsc/ldsc.py --h2 simphenos/res_meta.ldready \
	--ref-ld 1000g/simulated_ldscores \
	--w-ld 1000g/simulated_ldscores \
	--out simphenos/ldscm

