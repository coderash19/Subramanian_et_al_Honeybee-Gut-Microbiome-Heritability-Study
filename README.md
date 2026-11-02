### This Project contains bioinformatics and R statistical analysis pipelines for Subramanian et al. "Honey bee genetics shapes the abundance of core gut microbiome bacteria linked to health" (Submitted to Apidologie, 2026)

**Study Abstract:**
The gut microbiome of Western honey bees (Apis mellifera) is crucial for host health. While the environment is known to influence the composition of the honey bee microbiome, the role of genetics is not as well understood. Researchers have discovered differences in the microbiome between honey bee subspecies, but it is not currently known if genetics impact gut microbiome composition within colonies. Our study investigated the broad-sense heritability of the honey bee’s gut microbiome using patriline analysis. We used 16s metabarcoding to characterize the hindgut microbiome communities of 500 honey bee workers from 2 colonies. These workers were subsequently genotyped at several hypervariable microsatellite loci to assign workers into patrilineal groups. There was no evidence for heritability of the alpha diversity or beta diversity of the entire gut microbiome community, nor for the presence of most bacterial genera. However, we discovered that patriline identity significantly influenced the relative abundance of three ‘core’ bacteria genera (Gilliamella, Snodgrassella, and Bifidobacterium) that play key roles in modulating immune responses and detoxification of plant-derived compounds. Our findings suggest that, in addition to pharmacological interventions (i.e. Probiotics), selection for colonies with healthier gut microbiomes may be an avenue to improve the health of honey bees.

### Files for this Project Include: 

**FILE 1: Subramanian_et_al_BIOINFORMATICS_CODE.md**
The code for this part of the pipeline requires Linux/Bash

These Include:

**1. Patriline Assignment using 'COLONY' software**

**2. Microbiome Sequences Processing Pipeline**
- BASECALLING
- DEMULTIPLEXING
- ADAPTER AND PRIMER TRIMMING
- SIZE FILTERING
- INSTALLING AND PREPPING FOR QIIME2
- QIIME2 STEPS

**FILE 2: Subramanian_et_al_R_ANALYSIS.md**
The code for this part of the pipeline in this file requires R

These Include:
- CONVERTING QIIME QZA to CSV FILES FOR R
- MISTAG FILTERING
- DATA FILTERING AND WRANGLING
- RAREFACTION
- CONVERT TO RELATIVE ABUNDANCE and PRESENCE/ABSENCE TABLES
- GENERAL LINEAR MODELS AND PLOTS


