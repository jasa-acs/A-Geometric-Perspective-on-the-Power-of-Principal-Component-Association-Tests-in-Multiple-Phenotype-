# A Geometric Perspective on the Power of Principal Component Association Tests in Multiple Phenotype Studies

# Author Contributions Checklist Form

# Data

## Abstract
The data sets used in this paper are GWAS summary statistics (effect sizes, standard errors etc)
for eight phenotypes including Body Mass Index (BMI), Fasting Glucose (FG), Fasting Insulin
(FI), High-Density Lipoprotein cholesterol (HDL), Low-Density Lipoprotein cholesterol (LDL),
triglycerides (TG), Waist-hip-ratio (WHR), and Systolic Blood Pressure (SBP). We first merged
these single trait GWAS summary statistics data sets using the common 1,999,568 SNPs
shared by those single trait GWAS studies . We then performed our proposed PC based tests
using these univariate Z-scores.

## Availability
The data sets used in this paper are publicly availalable and the access links are provided in the
Description section.

## Description
The data sets used in this paper are publicly available and the links to download them are
provided below.

The BMI data set can be downloaded at
http://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz , GWAS Anthropometric 2015 BMI of European ancestry.

The reference paper is: Locke AE*, Kahali B*, Berndt SI*, Justice AE*, Pers TH*, Day FR,
Powell C, Vedantam S, Buchkovich ML, Yang J, Croteau-Chonka DC, Esko T et al. (2015)
Genetic studies of body mass index yield new insights for obesity biology. Nature 518, 197-206.

The WHR data set can be downloaded at
http://portals.broadinstitute.org/collaboration/giant/images/e/eb/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz

The reference paper is: Shungin D*, Winkler TW*, Croteau-Chonka DC*, Ferreira T*, Locke AE*,
Magi R*, Strawbridge R, Pers TH, Fischer K, Justice AE, Workalemahu T, Wu JM, et al. (2015)
New genetic loci link adipose and insulin biology to body fat distribution. Nature 518, 187-196.

The lipids (HDL, LDL, TG) data sets can be downloaded at
http://csg.sph.umich.edu/abecasis/public/lipids2013 .

The reference paper is: Willer CJ et al. Discovery and refinement of loci associated with lipid
levels. Nat. Genet. 2013. doi:10.1038/ng.2797

The MAGIC (FG and FI) data sets are at https://www.magicinvestigators.org/downloads/ in
the section of Glucose and insulin results accounting for BMI:

FG:
ftp://ftp.sanger.ac.uk/pub/magic/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt.gz

FI:
ftp://ftp.sanger.ac.uk/pub/magic/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt.gz

The reference paper is: Manning AK et al. A genome-wide approach accounting for body mass
index identifies genetic variants influencing fasting glycemic traits and insulin resistance.Nature
genetics 2012;44;6;659-69

The SBP data set is located at https://www.ncbi.nlm.nih.gov/projects/gap/cgibin/study.cgi?study_id=phs000585.v1.p1

The reference paper for SBP is:
Genetic variants in novel pathways influence blood pressure and cardiovascular disease risk.
Nature. 2011 Sep 11;478(7367):103-9.

# Code

## Abstract
The R code essentially merges the GWAS summary statistics data sets by the common
chromosome:position column and then apply the proposed PC based testing procedures to
each SNP. All the PC based tests are contained in the R package MPAT available at
https://content.sph.harvard.edu/xlin/dat/MPAT_1.0.tar.gz and will be available on CRAN
soon and github.

## Description
The main R package developed for the analysis is called MPAT publicly available at
https://content.sph.harvard.edu/xlin/dat/MPAT_1.0.tar.gz and the license is GPL (>=2).
The code will be avaiable on the github: https://github.com/zhonghualiu/GWAS .

## Optional Information 

The analysis can be done using modern personal computers with enought RAM to read in the
data (total size is about 560MB) or on high performance cluster.

# Instructions for Use

## Reproducibility

1. The subfolder named “Figure1_rejectionBoundary” contains the R code to draw
the rejection boundaries as shown in Figure 1.
2. The subfolder named “Figure2_Rotation” contains the R code to draw Figure 2.
3. The subfolder named “Figure3_Table5_Table6” contains the R code to draw
Figure 3 and generate the numbers for Table 5 and Table 6.
4. The subfolder named “Merge_analysis” contains the R code to merge all the raw
GWAS summary statistics for the eight phenotypes used in the paper. We merged
those GWAS summary statistics by using their chromosome: position column, and then
performed PC based tests to the Z-scores of each phenotype across all the SNPs in the
merged data. The output file is a csv file which contains the p-values for each PC based
method for each SNP.
5. The subfolder named “FigureS2-FigureS8_QQplotCode” contains the R code to
draw QQ plots and compute genomic control factors using the output csv file.
6. The subfolder named “FigureS1” contains the R code to draw Figure S1.
7. The subfolder named “SimulationSize” contains the R code to perform simulation
studies for the type I error rates, and the output was summarized in Table 1.
8. The subfolder named “SimulationPower” contains the R code to perform
simulation studies for the power comparisons, and the results were summarized in
Table 2 (simulation setup) and Table 3 (simulation results).
