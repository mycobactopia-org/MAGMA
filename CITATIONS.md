# TORCH-consortium/magma: Citations

## MAGMA pipeline

> Heupink TH, Verboven L, Sharma A, Rennie V, de Diego Fuertes M, Warren RM, Van Rie A.
> MAGMA: a comprehensive bioinformatics tool for genomic analysis of Mycobacterium tuberculosis.
> PLoS Comput Biol. 2023;19(12):e1011648. doi: 10.1371/journal.pcbi.1011648

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S.
> The nf-core framework for community-curated bioinformatics pipelines.
> Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C.
> Nextflow enables reproducible computational workflows.
> Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

### Read quality control

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

  > Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

  > Ewels P, Magnusson M, Lundin S, Käller M.
  > MultiQC: summarize analysis results for multiple tools and samples in a single report.
  > Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. PubMed PMID: 27312411.

### Read alignment

- [BWA](https://pubmed.ncbi.nlm.nih.gov/19451168/)

  > Li H, Durbin R.
  > Fast and accurate short read alignment with Burrows-Wheeler Aligner.
  > Bioinformatics. 2009 Jul 15;25(14):1754-60. doi: 10.1093/bioinformatics/btp324. PubMed PMID: 19451168.

- [SAMtools](https://pubmed.ncbi.nlm.nih.gov/19505943/)

  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup.
  > The Sequence Alignment/Map format and SAMtools.
  > Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. PubMed PMID: 19505943.

### Variant calling and genotyping

- [GATK](https://pubmed.ncbi.nlm.nih.gov/20644199/)

  > McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA.
  > The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data.
  > Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. PubMed PMID: 20644199.

  > Van der Auwera GA, Carneiro MO, Hartl C, Poplin R, Del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella KV, Altshuler D, Gabriel S, DePristo MA.
  > From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline.
  > Curr Protoc Bioinformatics. 2013;43:11.10.1-33. doi: 10.1002/0471250953.bi1110s43.

- [LoFreq](https://pubmed.ncbi.nlm.nih.gov/23066108/)

  > Wilm A, Aw PP, Bertrand D, Yeo GH, Ong SH, Wong CH, Khor CC, Petric R, Hibberd ML, Nagarajan N.
  > LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets.
  > Nucleic Acids Res. 2012 Dec;40(22):11189-201. doi: 10.1093/nar/gks918. PubMed PMID: 23066108.

### Structural variant detection

- [DELLY](https://pubmed.ncbi.nlm.nih.gov/22962449/)

  > Rausch T, Zichner T, Schlattl A, Stütz AM, Benes V, Korbel JO.
  > DELLY: structural variant discovery by integrated paired-end and split-read analysis.
  > Bioinformatics. 2012 Sep 15;28(18):i333-i339. doi: 10.1093/bioinformatics/bts378. PubMed PMID: 22962449.

- [ISMapper](https://pubmed.ncbi.nlm.nih.gov/26041500/)

  > Hawkey J, Edwards DJ, Dimopolous K, Hirani R, Phan MD, Tran N, Wahl C, Mulvey MR, Stinear TP, Dougan G, Holt KE.
  > ISMapper: identifying transposase insertion sites in bacterial genomes from short read sequence data.
  > BMC Genomics. 2015 Jun 4;16:329. doi: 10.1186/s12864-015-1660-1. PubMed PMID: 26041500.

### Variant annotation

- [SnpEff](https://pubmed.ncbi.nlm.nih.gov/22728672/)

  > Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM.
  > A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff.
  > Fly (Austin). 2012 Apr-Jun;6(2):80-92. doi: 10.4161/fly.19695. PubMed PMID: 22728672.

- [BCFtools](https://pubmed.ncbi.nlm.nih.gov/21903627/)

  > Li H.
  > A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data.
  > Bioinformatics. 2011 Nov 1;27(21):2987-93. doi: 10.1093/bioinformatics/btr509. PubMed PMID: 21903627.

### Drug resistance profiling

- [TBprofiler](https://pubmed.ncbi.nlm.nih.gov/33206998/)

  > Phelan JE, O'Sullivan DM, Machado D, Ramos J, Oppong YEA, Campino S, O'Grady J, McNerney R, Hibberd ML, Viveiros M, Huggett JF, Clark TG.
  > Integrating informatics tools and portable sequencing technology for rapid detection of resistance to anti-tuberculous drugs.
  > Genome Med. 2019 Feb 4;11(1):41. doi: 10.1186/s13073-019-0650-x. PubMed PMID: 33206998.

- [NTM-Profiler](https://pubmed.ncbi.nlm.nih.gov/35220982/)

  > Doyle RM, Poulain S, Walker TM, Fowler PW, Ashton P, Edwards A, Moore F, Smith EG, Satta G, Bhatt S, Crook DW, Wilson DJ, Harris SR.
  > Rapid and accurate species identification of Mycobacterium species from cultured isolates using whole genome sequencing.
  > J Clin Microbiol. 2022 Apr 20;60(4):e0200821. doi: 10.1128/jcm.02008-21. PubMed PMID: 35220982.

### Spoligotyping

- [SpoTyping](https://pubmed.ncbi.nlm.nih.gov/27234677/)

  > Xia E, Teo YY, Ong RT.
  > SpoTyping: fast and accurate in silico Mycobacterium spoligotyping from sequence reads.
  > Genome Med. 2016 May 27;8(1):19. doi: 10.1186/s13073-016-0270-7. PubMed PMID: 27234677.

### Phylogenetics

- [IQ-TREE 2](https://pubmed.ncbi.nlm.nih.gov/32011700/)

  > Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R.
  > IQ-TREE 2: New Models and Methods for Phylogenetic Inference.
  > Mol Biol Evol. 2020 May 1;37(5):1530-1534. doi: 10.1093/molbev/msaa015. PubMed PMID: 32011700.

- [snp-sites](https://pubmed.ncbi.nlm.nih.gov/28348851/)

  > Page AJ, Taylor B, Delaney AJ, Soares J, Seemann T, Keane JA, Harris SR.
  > SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments.
  > Microb Genom. 2016 Apr 29;2(4):e000056. doi: 10.1099/mgen.0.000056. PubMed PMID: 28348851.

- [snp-dists](https://github.com/tseemann/snp-dists)

  > Seemann T. snp-dists: Pairwise SNP distance matrix from a FASTA sequence alignment. GitHub. https://github.com/tseemann/snp-dists

- [ClusterPicker](https://pubmed.ncbi.nlm.nih.gov/24425371/)

  > Rambaut A, Lam TT, Max Carvalho L, Pybus OG.
  > Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen).
  > Virus Evol. 2016 Jan 1;2(1):vew007. doi: 10.1093/ve/vew007.

  > Rambaut A, Lam TT, Carvalho LM, Pybus OG. (ClusterPicker)
  > Cluster Picker: A Tool for Identifying Clusters in Phylogenetic Trees.
  > J Comput Biol. 2014 Jun;21(6):416-25. doi: 10.1089/cmb.2013.0220. PubMed PMID: 24425371.

### Compression and indexing

- [bgzip / htslib](https://pubmed.ncbi.nlm.nih.gov/33594379/)

  > Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H.
  > Twelve years of SAMtools and BCFtools.
  > Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PubMed PMID: 33594379.

## Reference datasets

- [EXIT-RIF lineage reference GVCF (LineagesAndOutgroupV2)](https://pubmed.ncbi.nlm.nih.gov/36319436/)

  > Heupink TH, Verboven L, Cuypers B, Van Rie A; TORCH Consortium.
  > Comprehensive, whole-genome sequencing based genotyping of clinical Mycobacterium tuberculosis (Mtb) using the EXIT-RIF reference collection.
  > Microb Genom. 2022 Nov;8(11). doi: 10.1099/mgen.0.000882. PubMed PMID: 36319436.

- [UVP exclusion loci / variant truth sets (Coll 2014, Coll 2018, Walker 2015, Napier 2020, Benavente 2015, Zeng 2018)](https://pubmed.ncbi.nlm.nih.gov/28814674/)

  > Coll F, McNerney R, Preston MD, Guerra-Assunção JA, Warry A, Hill-Cawthra G, Mallard K, Nair M, Miranda A, Alves A, Perdigão J, Viveiros M, Portugal I, Hasan Z, Hasan R, Glynn JR, Martin N, Pain A, Parkhill J, McNerney R, Clark TG.
  > Rapid determination of anti-tuberculosis drug resistance from whole-genome sequences.
  > Genome Med. 2015 Sep 1;7(1):51. doi: 10.1186/s13073-015-0164-0.

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team.
  > Bioconda: sustainable and comprehensive software distribution for the life sciences.
  > Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y.
  > BioContainers: an open-source and community-driven framework for software standardization.
  > Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW.
  > Singularity: Scientific containers for mobility of compute.
  > PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. PubMed PMID: 28494014.
