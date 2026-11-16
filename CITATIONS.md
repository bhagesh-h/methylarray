# methylarray: Citations

This file is the single reference list for the pipeline. The README, the documentation and the reports point here instead of repeating references. The first sections list the software, methods and data the pipeline runs or reads; the last sections list work that informed its design but is not run by it. Entries without a DOI link have no DOI or are pending verification of the identifier.

## Pipeline

- methylarray

  > A Zenodo DOI will be minted at the first release. Until then, cite the GitHub repository and the exact version or commit you ran.

## Frameworks

- [Nextflow](https://doi.org/10.1038/nbt.3820)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017;35(4):316-319.

- [nf-core](https://doi.org/10.1038/s41587-020-0439-x)

  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020;38(3):276-278.

  The pipeline is built from the nf-core template and reuses nf-core code and infrastructure under the MIT licence.

- [nf-schema](https://github.com/nextflow-io/nf-schema)

  > Parameter and samplesheet validation plugin for Nextflow.

- [nf-test](https://www.nf-test.com)

  > Forer L, Schönherr S. Improving the reliability and quality of Nextflow pipelines with nf-test. bioRxiv. 2024.

- [MultiQC](https://doi.org/10.1093/bioinformatics/btw354)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048.

- [Bioconductor](https://doi.org/10.1038/nmeth.3252)

  > Huber W, Carey VJ, Gentleman R, et al. Orchestrating high-throughput genomic analysis with Bioconductor. Nat Methods. 2015;12(2):115-121.

- [Quarto](https://quarto.org)

  > Open-source scientific and technical publishing system. Renders the HTML analysis report from `assets/report/methylarray_report.qmd`.

## Software packaging and containers

- [Bioconda](https://doi.org/10.1038/s41592-018-0046-7)

  > Grüning B, Dale R, Sjödin A, et al. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018;15(7):475-476.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel D. Docker: lightweight Linux containers for consistent development and deployment. Linux Journal. 2014;2014(239):2.

- [Singularity](https://doi.org/10.1371/journal.pone.0177459)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: scientific containers for mobility of compute. PLoS One. 2017;12(5):e0177459.

- [Seqera Containers and Wave](https://seqera.io/containers/)

  > Container images built on demand from Conda environment files.

- [nf-metro](https://github.com/seqeralabs/nf-metro)

  > Metro-map rendering of the workflow diagram (`docs/images/metro_map.svg`) from `assets/metro_map.mmd`.

## Array annotation and probe characterisation

- [EPIC v2 characterisation and enhanced manifest](https://doi.org/10.1186/s12864-024-10027-5)

  > Peters TJ, Meyer B, Ryan L, et al. Characterisation and reproducibility of the HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling. BMC Genomics. 2024;25:251.

- [Probe masking annotation](https://doi.org/10.1093/nar/gkw967)

  > Zhou W, Laird PW, Shen H. Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Res. 2017;45(4):e22.

- [EPICv2manifest](https://doi.org/10.18129/B9.bioc.EPICv2manifest)

  > Peters T. EPICv2manifest: Illumina Infinium MethylationEPIC v2.0 extended manifest from Peters et al. 2024. Bioconductor package. The pipeline reads AnnotationHub resource AH116484.

- [InfiniumAnnotation](https://github.com/zhou-lab/InfiniumAnnotationV1)

  > Zhou lab probe annotation and masks for Infinium arrays. The pipeline reads the EPIC v2 hg38 mask table (`EPICv2.hg38.mask.tsv.gz`) at a fixed commit.

- IlluminaHumanMethylationEPICv2manifest and IlluminaHumanMethylationEPICv2anno.20a1.hg38

  > Bioconductor annotation packages for the EPIC v2 array (hg38).

## Preprocessing and quality control

- [SeSAMe](https://doi.org/10.1093/nar/gky691)

  > Zhou W, Triche TJ Jr, Laird PW, Shen H. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. Nucleic Acids Res. 2018;46(20):e123.

- [noob background correction](https://doi.org/10.1093/nar/gkt090)

  > Triche TJ Jr, Weisenberger DJ, Van Den Berg D, Laird PW, Siegmund KD. Low-level processing of Illumina Infinium DNA Methylation BeadArrays. Nucleic Acids Res. 2013;41(7):e90.

- [minfi](https://doi.org/10.1093/bioinformatics/btu049)

  > Aryee MJ, Jaffe AE, Corrada-Bravo H, et al. Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. Bioinformatics. 2014;30(10):1363-1369.

- [Functional normalisation](https://doi.org/10.1186/s13059-014-0503-2)

  > Fortin JP, Labbe A, Lemire M, et al. Functional normalization of 450k methylation array data improves replication in large cancer studies. Genome Biol. 2014;15:503.

- [SWAN normalisation](https://doi.org/10.1186/gb-2012-13-6-r44)

  > Maksimovic J, Gordon L, Oshlack A. SWAN: Subset-quantile within array normalization for Illumina Infinium HumanMethylation450 BeadChips. Genome Biol. 2012;13(6):R44.

- [Stratified quantile normalisation](https://doi.org/10.2217/epi.12.21)

  > Touleimat N, Tost J. Complete pipeline for Infinium Human Methylation 450K BeadChip data processing using subset quantile normalization for accurate DNA methylation estimation. Epigenomics. 2012;4(3):325-341.

- [minfi for EPIC](https://doi.org/10.1093/bioinformatics/btw691)

  > Fortin JP, Triche TJ Jr, Hansen KD. Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array with minfi. Bioinformatics. 2017;33(4):558-560.

- [ENmix](https://doi.org/10.1093/nar/gkv907)

  > Xu Z, Niu L, Li L, Taylor JA. ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip. Nucleic Acids Res. 2016;44(3):e20.

  > Xu Z, Niu L, Taylor JA. The ENmix DNA methylation analysis pipeline for Illumina BeadChip and comparisons with seven other preprocessing pipelines. Clin Epigenetics. 2021;13:216.

- [RCP probe-type correction](https://doi.org/10.1093/bioinformatics/btw285)

  > Niu L, Xu Z, Taylor JA. RCP: a novel probe design bias correction method for Illumina Methylation BeadChip. Bioinformatics. 2016;32(17):2659-2663.

- [RELIC dye-bias correction](https://doi.org/10.1186/s12864-016-3426-3)

  > Xu Z, Langie SA, De Boever P, Taylor JA, Niu L. RELIC: a novel dye-bias correction method for Illumina Methylation BeadChip. BMC Genomics. 2017;18:4.

- [ewastools](https://doi.org/10.1186/s13148-018-0504-1)

  > Heiss JA, Just AC. Identifying mislabeled and contaminated DNA methylation microarray data: an extended quality control toolset with examples from GEO. Clin Epigenetics. 2018;10:73.

- [omicsPrint (identity-by-state plot of array pairs)](https://doi.org/10.1093/bioinformatics/bty062)

  > van Iterson M, Cats D, Hop P, BIOS Consortium, Heijmans BT. omicsPrint: detection of data linkage errors in multiple omics studies. Bioinformatics. 2018;34(12):2142-2143.

- [Pairwise relationship likelihoods](https://doi.org/10.1111/j.1469-1809.1975.tb00120.x)

  > Thompson EA. The estimation of pairwise relationships. Ann Hum Genet. 1975;39(2):173-188.

- [illuminaio](https://doi.org/10.12688/f1000research.2-264.v1)

  > Smith ML, Baggerly KA, Bengtsson H, Ritchie ME, Hansen KD. illuminaio: An open source IDAT parsing tool for Illumina microarrays. F1000Research. 2013;2:264.

- [EpiAnceR+ (genetic ancestry components from methylation arrays)](https://doi.org/10.1186/s13072-025-00627-0)

  > Höffler KD, Katrinli S, Halvorsen MW, et al. Optimizing genetic ancestry adjustment in DNA methylation studies: a comparative analysis of approaches. Epigenetics Chromatin. 2025;18(1):69.

  The pipeline implements the published method; it does not use the EpiAnceR scripts.

- [GoDMC methylation QTLs](https://doi.org/10.1038/s41588-021-00923-x)

  > Min JL, Hemani G, Hannon E, et al. Genomic and phenotypic insights from an atlas of genetic effects on DNA methylation. Nat Genet. 2021;53(9):1311-1321.

  Source of the meQTL table described in the usage documentation (`--meqtl`); the pipeline does not ship it.

- [Genotype-robust clocks](https://doi.org/10.1038/s42003-025-08893-0)

  > Meeks GL, Scelza B, Asnake HM, et al. Common DNA sequence variation influences epigenetic aging in African populations. Commun Biol. 2025;8(1).

- [variancePartition](https://doi.org/10.1186/s12859-016-1323-z)

  > Hoffman GE, Schadt EE. variancePartition: interpreting drivers of variation in complex gene expression studies. BMC Bioinformatics. 2016;17:483.

- [Control-probe principal components](https://doi.org/10.1186/s13059-015-0600-x)

  > Lehne B, Drong AW, Loh M, et al. A coherent approach for analysis of the Illumina HumanMethylation450 BeadChip improves data quality and performance in epigenome-wide association studies. Genome Biol. 2015;16:37.

- [isva (random matrix theory estimate of the number of components)](https://doi.org/10.1093/bioinformatics/btr171)

  > Teschendorff AE, Zhuang J, Widschwendter M. Independent surrogate variable analysis to deconvolve confounding factors in large-scale microarray profiling studies. Bioinformatics. 2011;27(11):1496-1505.

- SummarizedExperiment and HDF5Array

  > Bioconductor infrastructure packages for on-disk assay storage.

- sesameData and ExperimentHub

  > Bioconductor packages that distribute the SeSAMe probe manifests, masks and platform signatures. The `sesame` and `combined` routes download them once per run (or reuse `--sesame_cache`) and records each resource's ExperimentHub ID and MD5 checksum.

## Cell-type composition

- [EpiDISH](https://doi.org/10.1186/s12859-017-1511-5)

  > Teschendorff AE, Breeze CE, Zheng SC, Beck S. A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies. BMC Bioinformatics. 2017;18:105.

- [12-cell-type blood reference (cent12CT)](https://doi.org/10.1186/s13073-023-01211-5)

  > Luo Q, Dwaraka VB, Chen Q, et al. A meta-analysis of immune-cell fractions at high resolution reveals novel associations with common phenotypes and health outcomes. Genome Med. 2023;15:59.

## DNA methylation biomarkers

- [Horvath pan-tissue clock](https://doi.org/10.1186/gb-2013-14-10-r115)

  > Horvath S. DNA methylation age of human tissues and cell types. Genome Biol. 2013;14:R115.

- [Hannum clock](https://doi.org/10.1016/j.molcel.2012.10.016)

  > Hannum G, Guinney J, Zhao L, et al. Genome-wide methylation profiles reveal quantitative views of human aging rates. Mol Cell. 2013;49(2):359-367.

- [Skin and blood clock](https://doi.org/10.18632/aging.101508)

  > Horvath S, Oshima J, Martin GM, et al. Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome and ex vivo studies. Aging. 2018;10(7):1758-1775.

- [PhenoAge](https://doi.org/10.18632/aging.101414)

  > Levine ME, Lu AT, Quach A, et al. An epigenetic biomarker of aging for lifespan and healthspan. Aging. 2018;10(4):573-591.

- [GrimAge (basis of PC-GrimAge)](https://doi.org/10.18632/aging.101684)

  > Lu AT, Quach A, Wilson JG, et al. DNA methylation GrimAge strongly predicts lifespan and healthspan. Aging. 2019;11(2):303-327.

- [PC clocks](https://doi.org/10.1038/s43587-022-00248-2)

  > Higgins-Chen AT, Thrush KL, Wang Y, et al. A computational solution for bolstering reliability of epigenetic clocks: implications for clinical trials and longitudinal tracking. Nat Aging. 2022;2:644-661.

- [DunedinPACE](https://doi.org/10.7554/eLife.73420)

  > Belsky DW, Caspi A, Corcoran DL, et al. DunedinPACE, a DNA methylation biomarker of the pace of aging. eLife. 2022;11:e73420.

- [dnaMethyAge](https://doi.org/10.1007/s11357-023-00871-w) ([GitHub](https://github.com/yiluyucheng/dnaMethyAge))

  > Wang Y, Grant OA, Zhai X, McDonald-Maier KD, Schalkwyk LC. Insights into ageing rates comparison across tissues from recalibrating cerebellum DNA methylation clock. GeroScience. 2024;46:39-56.

  Implementation used for the published clock coefficients (pinned to a commit). It also provides the clocks below; cAge (BernabeuE2023c) is not in the default `--clocks` list.

- Zhang 2019 clock

  > Zhang Q, Vallerga CL, Walker RM, et al. Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing. Genome Med. 2019;11:54.

- cAge

  > Bernabeu E, McCartney DL, Gadd DA, et al. Refining epigenetic prediction of chronological and biological age. Genome Med. 2023;15:12.

- DNAmTL

  > Lu AT, Seeboth A, Tsai PC, et al. DNA methylation-based estimator of telomere length. Aging. 2019;11(16):5895-5923.

- Mortality risk score

  > Zhang Y, Wilson R, Heiss J, et al. DNA methylation signatures in peripheral blood strongly predict all-cause mortality. Nat Commun. 2017;8:14617.

- epiTOC2

  > Teschendorff AE. A comparison of epigenetic mitotic-like clocks for cancer risk prediction. Genome Med. 2020;12:56.

- DNAmIC

  > Fuentealba M, Rouch L, Guyonnet S, et al. A blood-based epigenetic clock for intrinsic capacity predicts mortality and is associated with clinical, immunological and lifestyle factors. Nat Aging. 2025.

- [Smoking, alcohol and BMI predictors](https://doi.org/10.1186/s13059-018-1514-1)

  > McCartney DL, Hillary RF, Stevenson AJ, et al. Epigenetic prediction of complex traits and death. Genome Biol. 2018;19:136.

- [DNA methylation IL-6 score](https://doi.org/10.1093/gerona/glab046)

  > Stevenson AJ, Gadd DA, Hillary RF, et al. Creating and validating a DNA methylation-based proxy for interleukin-6. J Gerontol A Biol Sci Med Sci. 2021;76(12):2284-2292.

  The smoking, alcohol, BMI and IL-6 weights in `assets/dnam_scores` are taken from the supplementary tables of McCartney et al. (2018) and Stevenson et al. (2021), both published under CC BY 4.0; `assets/dnam_scores/sources.tsv` names each table.

- [Clocks in genetically admixed individuals](https://doi.org/10.7554/eLife.105343)

  > Cruz-Gonzalez S, Okpala O, Gu E, et al. Methylation clocks fail to generalize across genetically admixed individuals. eLife. 2026;14.

  Cited by the report when it relates the accelerations to the genetic ancestry components.

## Statistical analysis

- [limma](https://doi.org/10.1093/nar/gkv007)

  > Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47.

- [bacon](https://doi.org/10.1186/s13059-016-1131-9)

  > van Iterson M, van Zwet EW, BIOS Consortium, Heijmans BT. Controlling bias and inflation in epigenome- and transcriptome-wide association studies using the empirical null distribution. Genome Biol. 2017;18:19.

- [Beta values and M-values](https://doi.org/10.1186/1471-2105-11-587)

  > Du P, Zhang X, Huang CC, et al. Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. BMC Bioinformatics. 2010;11:587.

- [EPIC significance threshold](https://doi.org/10.1186/s12864-019-5761-7)

  > Mansell G, Gorrie-Stone TJ, Bao Y, et al. Guidance for DNA methylation studies: statistical insights from the Illumina EPIC array. BMC Genomics. 2019;20:366.

- [Effective number of tests](https://doi.org/10.1002/gepi.20408)

  > Galwey NW. A new measure of the effective number of tests, a practical tool for comparing families of non-independent significance tests. Genet Epidemiol. 2009;33(7):559-568.

## Differentially methylated regions and enrichment

- [DMRcate](https://doi.org/10.1186/1756-8935-8-6)

  > Peters TJ, Buckley MJ, Statham AL, et al. De novo identification of differentially methylated regions in the human genome. Epigenetics Chromatin. 2015;8:6.

- [dmrff](https://doi.org/10.1101/508556)

  > Suderman M, Staley JR, French R, Arathimos R, Simpkin A, Tilling K. dmrff: identifying differentially methylated regions efficiently with power and control. bioRxiv. 2018.

- [KnowYourCG](https://doi.org/10.1126/sciadv.adw3027)

  > Goldberg DC, Fu H, Atkins D, Moyer E, Lee CN, Deng Y, Zhou W. KnowYourCG: facilitating base-level sparse methylome interpretation. Sci Adv. 2025;11(43).

  The EPIC v2 knowledgebases are ExperimentHub records EH8068 to EH8075 (sesameData).

- [missMethyl](https://doi.org/10.1093/bioinformatics/btv560)

  > Phipson B, Maksimovic J, Oshlack A. missMethyl: an R package for analyzing data from Illumina's HumanMethylation450 platform. Bioinformatics. 2016;32(2):286-288.

- [Gene set testing for methylation data](https://doi.org/10.1186/s13059-021-02388-x)

  > Maksimovic J, Oshlack A, Phipson B. Gene set enrichment analysis for genome-wide DNA methylation data. Genome Biol. 2021;22:173.

- [MSigDB Hallmark gene sets](https://doi.org/10.1016/j.cels.2015.12.004)

  > Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database hallmark gene set collection. Cell Syst. 2015;1(6):417-425.

- [Reactome](https://doi.org/10.1093/nar/gkad1025)

  > Milacic M, Beavers D, Conley P, et al. The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res. 2024;52(D1):D672-D678.

- [KEGG](https://doi.org/10.1093/nar/28.1.27)

  > Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 2000;28(1):27-30.

- [Gene Ontology](https://doi.org/10.1093/genetics/iyad031)

  > The Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023. Genetics. 2023;224(1):iyad031.

## Public data

- GEO [GSE286313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286313): matched EPIC v1 and v2 venous blood (Zhuang et al. 2025). The EPIC v2 arrays are the data of the `test` and `test_full` profiles and of the module tests.

- GEO [GSE240469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240469): EPIC v2 cell line replicates (Peters et al. 2024).

## Related methods and resources not run by the pipeline

These informed the design or are candidates for later versions. The pipeline does not run them.

- [EPIC v2 evaluation](https://doi.org/10.1186/s43682-023-00021-5)

  > Kaur D, Lee SM, Goldberg D, et al. Comprehensive evaluation of the Infinium human MethylationEPIC v2 BeadChip. Epigenetics Commun. 2023;3:6.

- EPIC v2 validation

  > Noguera-Castells A, García-Prieto CA, Álvarez-Errico D, Esteller M. Validation of the new EPIC DNA methylation microarray (900K EPIC v2) for high-throughput profiling of the human DNA methylome. Epigenetics. 2023;18(1):2185742.

- [EPIC v1 probe evaluation](https://doi.org/10.1186/s13059-016-1066-1)

  > Pidsley R, Zotenko E, Peters TJ, et al. Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling. Genome Biol. 2016;17:208.

- [EPIC v1 polymorphic and off-target probes](https://doi.org/10.1016/j.gdata.2016.05.012)

  > McCartney DL, Walker RM, Morris SW, McIntosh AM, Porteous DJ, Evans KL. Identification of polymorphic and off-target probe binding sites on the Illumina Infinium MethylationEPIC BeadChip. Genom Data. 2016;9:22-24.

- [BMIQ](https://doi.org/10.1093/bioinformatics/bts680)

  > Teschendorff AE, Marabita F, Lechner M, et al. A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450 k DNA methylation data. Bioinformatics. 2013;29(2):189-196.

- [BioContainers](https://doi.org/10.1093/bioinformatics/btx192)

  > da Veiga Leprevost F, Grüning BA, Alves Aflitos S, et al. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017;33(16):2580-2582.

- [Enhanced blood deconvolution reference](https://doi.org/10.1038/s41467-021-27864-7)

  > Salas LA, Zhang Z, Koestler DC, et al. Enhanced cell deconvolution of peripheral blood using DNA methylation for high-resolution immune profiling. Nat Commun. 2022;13:761.

- [FlowSorted.Blood.EPIC and IDOL](https://doi.org/10.1186/s13059-018-1448-7)

  > Salas LA, Koestler DC, Butler RA, et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol. 2018;19:64.

- [Houseman constrained projection](https://doi.org/10.1186/1471-2105-13-86)

  > Houseman EA, Accomando WP, Koestler DC, et al. DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics. 2012;13:86.

- [methylCIPHER](https://doi.org/10.1101/2022.07.13.499978)

  > Thrush KL, Higgins-Chen AT, Liu Z, Levine ME. R methylCIPHER: a methylation clock investigational package for hypothesis-driven evaluation and research. bioRxiv. 2022.

- pyaging

  > de Lima Camillo LP. pyaging: a Python-based compendium of GPU-optimized aging clocks. Bioinformatics. 2024.

- [DNA methylation CRP score](https://doi.org/10.1038/s41467-022-29792-6)

  > Wielscher M, Mandaviya PR, Kuehnel B, et al. DNA methylation signature of chronic low-grade inflammation and its role in cardio-respiratory diseases. Nat Commun. 2022;13:2408.

- [EpiScores](https://doi.org/10.7554/eLife.71802)

  > Gadd DA, Hillary RF, McCartney DL, et al. Epigenetic scores for the circulating proteome as tools for disease prediction. eLife. 2022;11:e71802.

- [variancePartition](https://doi.org/10.1186/s12859-016-1323-z)

  > Hoffman GE, Schadt EE. variancePartition: interpreting drivers of variation in complex gene expression studies. BMC Bioinformatics. 2016;17:483.

- [sva](https://doi.org/10.1093/bioinformatics/bts034)

  > Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD. The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics. 2012;28(6):882-883.

- [ComBat](https://doi.org/10.1093/biostatistics/kxj037)

  > Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007;8(1):118-127.

- [450K significance threshold](https://doi.org/10.1002/gepi.22086)

  > Saffari A, Silver MJ, Zavattari P, et al. Estimation of a significance threshold for epigenome-wide association studies. Genet Epidemiol. 2018;42(1):20-33.

- [Batch correction caveats](https://doi.org/10.3389/fgene.2018.00083)

  > Price EM, Robinson WP. Adjusting for batch effects in DNA methylation microarray data, a lesson learned. Front Genet. 2018;9:83.

- [ComBat simulation study](https://doi.org/10.1186/s12859-020-03559-6)

  > Zindler T, Frieling H, Neyazi A, Bleich S, Friedel E. Simulating ComBat: how batch correction can lead to the systematic introduction of false positive results in DNA methylation microarray studies. BMC Bioinformatics. 2020;21:271.

- KnowYourCG

  > CpG-centric enrichment testing distributed with SeSAMe (see the SeSAMe citation).

## Reference studies for validating results

- [EWAS of body mass index](https://doi.org/10.1038/nature20784)

  > Wahl S, Drong A, Lehne B, et al. Epigenome-wide association study of body mass index, and the adverse outcomes of adiposity. Nature. 2017;541:81-86.

- The EWAS Catalog

  > Battram T, Yousefi P, Crawford G, et al. The EWAS Catalog: a database of epigenome-wide association studies. Wellcome Open Res. 2022;7:41.

## Background reviews that informed the design

- [Analysing and interpreting DNA methylation data](https://doi.org/10.1038/nrg3273)

  > Bock C. Nat Rev Genet. 2012;13(10):705-719.

- [A practical guide to the measurement and analysis of DNA methylation](https://doi.org/10.1165/rcmb.2019-0150TR)

  > Singer BD. Am J Respir Cell Mol Biol. 2019;61(4):417-428.

- [Methylation data analysis and interpretation](https://doi.org/10.1146/annurev-biodatasci-120924-091033)

  > Zhu Y, Mao W, Hosseini R, Chikina M. Annu Rev Biomed Data Sci. 2025.

- [Differences between EPIC v2 and v1 in DNA methylation based tools](https://doi.org/10.26508/lsa.202403155)

  > Zhuang BC, Jude MS, Konwar C, et al. Accounting for differences between Infinium MethylationEPIC v2 and v1 in DNA methylation-based tools. Life Sci Alliance. 2025;8(9):e202403155.
