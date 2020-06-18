# Version1
FPP-FM: a general framework for functional mapping of physiological big data from phenomics assay


This package is developed to detect quantitative trait loci (QTL) for multi-stage physiological big data from phenomics assay based on the Functional mapping [1], these phenotype data is measured by high-throughput phenotyping techniques. This manual describes some instructions on how to perform the tasks of QTL detection by this package. The briefing of this guide is as follow: 

Section 1: Data format

Section 2: Examples

#1. Data format

#This package works on two files include genotype and phenotype. The appropriate formatting for these two files is described below.

#Genotype file

#ID	SNP1	SNP2	SNP3	SNP4	SNP5	SNP6	SNP7	...

#1	1	1	1	0	0	0	0	...

#2	0	1	1	0	1	1	0	...

#3	0	0	1	1	0	1	0	...

#4	0	1	1	1	1	0	1	...

#5	1	0	0	1	1	1	0	...

#...								

#Two genotypes (QQ=1, qq=0).

#Phenotype file

#ID	Point1	Point2	Point3	Point4	Point5	Point6	Point7	...

#1	0.125	0.289	0.308	0.232	0.383	0.237	0.208	...

#2	0.184	0.397	0.520	0.382	0.454	0.372	0.213	...

#3	0.151	0.249	0.413	0.361	0.420	0.367	0.249	...

#4	0.144	0.371	0.552	0.418	0.443	0.356	0.301	...

#5	0.150	0.254	0.284	0.265	0.367	0.262	0.292	...

#...						

#The phenotype file contains numerical values.

#2. Example

#set working directory

setwd("~/code/")

#########Load mvtnorm package and functions###############

library(mvtnorm) #Available on CRAN

source("FP_sup.R") 

######## Load data for genotype file and phenotype file ###################

dat <- data.load(pheno="pheno.csv",marker="genotype.csv",time=1:31)

#fw.dat.load function has three arguments, including # 1)pheno indicates phenotype file; 2) marker indicates genotype file; 3)time indicates time point of the actual measurement.

#dat list includes the following:

#1) Phenotype table: phenotype.

#2) Phenotype table: genotype.

#3) Vector (time point or other indicator variable): time.

######## Null hypothesis#########################

H0 <- mle_curve(pheno=dat$phenotype,times=dat$time)

#H0 is a list include the initial parameters of covariance matrix and Legendre's polynomiials.

######## Genome scan ############################

ret <- mle_H1(dat,times=dat$time)

#ret matrix includes LR , L1and parameters of covariance matrix and equations. Each row of the ret matrix is the result of hypothesis test for each SNP, similar to the initial parameters. The second column is the LR value of all SNP, which can be map Manhattan plot. 

#[1] Ma, C. X., Casella, G., & Wu, R. (2002). Functional mapping of quantitative trait loci underlying the character process: a theoretical framework. Genetics, 161(4), 1751-1762.



