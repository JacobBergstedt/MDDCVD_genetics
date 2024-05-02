singularity exec --home $PWD $SIF/ldsc.sif python /tools/ldsc/munge_sumstats.py \
 --sumstats "Data/Sumstats/MDDCVD_Q0.05_latent_factor_sumstats" \
 --out "./Data/Munged_sumstats/MDDCVD_Q0.05_latent_factor_sumstats_munged"\
 --chunksize 500000 \
 --signed-sumstats Z,0 \
 --merge-alleles /REF/ldsc/w_hm3.snplist

singularity exec --home $PWD $SIF/ldsc.sif python /tools/ldsc/ldsc.py \
--h2 "Data/Munged_sumstats/MDDCVD_Q0.05_latent_factor_sumstats_munged.sumstats" \
--ref-ld-chr /REF/ldsc/eur_w_ld_chr/ \
--w-ld-chr /REF/ldsc/eur_w_ld_chr/ \
--out "Results/gSEM_joelle/MDDCVD_Q0.05_h2"

singularity exec --home $PWD $COMORMENT/ldsc/containers/ldsc.sif python /tools/ldsc/ldsc.py \
--h2 "Data/Munged_sumstats/CVD_MDD_GWAS_V2.5_formatted.txt_munged.sumstats.gz" \
--ref-ld-chr Data/REF/eur_w_ld_chr/ \
--w-ld-chr Data/REF/eur_w_ld_chr/ \
--out "Results/gSEM_joelleCVD_MDD_GWAS_V2.5_h2"

singularity exec --home $PWD $COMORMENT/ldsc/containers/ldsc.sif python /tools/ldsc/ldsc.py \
--h2 "Data/Munged_sumstats/CVD_MDD_GWAS_V2.5_Q0.05_formatted.txt_munged.sumstats" \
--ref-ld-chr Data/REF/eur_w_ld_chr/ \
--w-ld-chr Data/REF/eur_w_ld_chr/ \
--out "Results/gSEM_joelleCVD_MDD_GWAS_V2.5_Q0.05_h2"


### Update

# Before filtering

singularity exec --home $PWD $COMORMENT/ldsc/containers/ldsc.sif python /tools/ldsc/ldsc.py \
--h2 "Data/Munged_sumstats/MDDCVD_Als_2023_latent_factor_Q0.05_sumstats_munged.sumstats" \
--ref-ld-chr Data/REF/eur_w_ld_chr/ \
--w-ld-chr Data/REF/eur_w_ld_chr/ \
--out "Results/gSEM_joelle/CVD_MDD_GWAS_Als_Q0.05_h2"

singularity exec --home $PWD $COMORMENT/ldsc/containers/ldsc.sif python /tools/ldsc/ldsc.py \
--h2 "Data/Munged_sumstats/MDDCVD_Als_2023_latent_factor_sumstats_munged.sumstats" \
--ref-ld-chr Data/REF/eur_w_ld_chr/ \
--w-ld-chr Data/REF/eur_w_ld_chr/ \
--out "Results/gSEM_joelle/CVD_MDD_GWAS_Als_h2"
