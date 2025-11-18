# MACHINE-paper

This repository contains the analysis code and results supporting the findings of the MACHINE manuscript (https://www.researchsquare.com/article/rs-7737326/v1).
The provided materials include:

- Scripts for data processing and analysis.

- All scripts and minimal datasets used to reproduce the manuscript's results and visualizations, which are documented in `data_info.xlsx`.

For detailed instructions on the analyses, refer to the sections below.


## **simulation**

This directory contains scripts for generating LD matrices, simulating summary statistics, performing fine-mapping analyses, visualization, and summarized results of simulation studies. Below is a detailed description of each script.

### **1. LD matrix generation**

- **simulation_data_LD.R**

    Select 200 LD blocks and generate LD matrices for EUR and EAS populations from UK BioBank (UKBB) genotype data.

- **simulation_data_LD_1kg.R**

    Compute LD matrices for EUR and EAS populations using 1000 Genomes (1kG) genotype data.


### **2. Simulated summary statistics**

- **simulation_data.R**

    Generate z-scores from the UKBB LD matrices.


### **3. Functional annotation integration**

- **do_gLDSC.R**

    Run g-LDSC for each simulation setting.

- **gLDSC_var.R**

    Output per-variant heritability estimated by g-LDSC.

- **polyfun_input.R**

    Prepares input data for PolyFun.

- **do_polyfun.sh**

    Execute PolyFun for each simulation setting (under polyfun environment).

- **polyfun_var.R**

    Output per-variant heritability estimated by PolyFun.


### **4. Fine mapping with UKBB LD matrices**

- **simulation_[method].R**

    Perform fine mapping using UKBB LD matrices with the following methods: 
    
    + MACHINE-anno (MACHINE + g-LDSC / MACHINE + PolyFun)
    + MACHINE
    + MESuSiE
    + SuSiEx
    + XMAP
    + MultiSuSiE
    + h2D2-anno (h2-D2 + g-LDSC / h2-D2 + PolyFun)
    + h2D2 (h2-D2)
    + SuSiE
    + CARMA


### **5. Fine mapping with UKBB LD matrices**

- **simulation\_[method]\_1kg.R**

    Perform fine mapping using 1kG LD matrices with the following methods:
    
    + MACHINE-anno (MACHINE + g-LDSC / MACHINE + PolyFun)
    + MACHINE
    + h2D2-anno (h2-D2 + g-LDSC / h2-D2 + PolyFun)
    + h2D2 (h2-D2)
    + RSparsePro
    + CARMA



## **real_data**

The `glgc` directory contains scripts for processing GWAS summary statistics, performing fine-mapping analyses, and summarized fine-mapping results for four lipid traits using data from UKBB and the Global Lipids Genetics Consortium (GLGC).

The `scz2022` directory contains scripts for processing GWAS summary statistics, performing fine-mapping analyses, and summarized fine-mapping results for schizophrenia (SCZ) using data from Psychiatric Genomics Consortium (PGC).


### 1. **Data preprocessing**

- **data_preprocess.R**

    Preprocess raw GWAS summary statistics.

- **glgc/data_merge_[db].R**

    Merge multi-ancestry summary statistics for UKBB and GLGC datasets.

- **glgc/min_p_sub.R**

    Filter loci containing at least one variant with a marginal association $P < 10^{-5}$ in any ancestry.

- **scz2022/data_merge.R**

    Merge multi-ancestry summary statistics. Filter loci containing at least one variant with a marginal association $P < 10^{-5}$ in any ancestry.

- **LD_overlap.R**

    Preprocess LD matrices from the UKBB reference panel.


### 2. **Functional annotation integration**

- **do_gLDSC.R**

    Perform g-LDSC for each trait and each ancestry.

- **gLDSC_var.R**

    Output per-variant heritability estimated by g-LDSC for each trait and each ancestry.


### 3. **Fine mapping**

- **do_[method].R**

    Perform fine mapping using the following methods:

    + MACHINE-gLDSC (MACHINE + g-LDSC)
    + MACHINE
    + MESuSiE
    + h2D2-gLDSC (h2-D2 + g-LDSC), 
    + h2D2 (h2-D2)
    + SuSiE
    + RSparsePro
    + CARMA


### 4. **External functional annotation preprocessing**

- **dbSNP151_preprocess.R**

    Preprocess gene-based annotations from dbSNP151.

- **GTEx_eQTL_preprocess.R**

    Preprocess GTEx eQTL data.


### 5. **Results**

- **results/**

    Directory containing summarized fine-mapping results.


## **gLDSC**

This directory contains scripts for performing g-LDSC. For more information, please refer to https://github.com/xzw20046/gldsc.