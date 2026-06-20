# References — methods and tools used in this pipeline

The tool citations below were retrieved from **PubMed** (via the bio-research PubMed
connector). Each entry lists the paper that introduces the method and the role it plays
in this pipeline. DOIs are linked. A few tools (STAR-Fusion, Picard, Nextflow, nf-core,
AMP/ASCO/CAP) are cited from their primary publications/software pages where a PubMed
record was not the canonical source.

> Attribution: tool papers in the "PubMed-sourced" list were obtained from PubMed; please
> keep the DOI links when reusing this list.

## Core methods (PubMed-sourced)

1. **STAR — spliced alignment.** Dobin A, Davis CA, Schlesinger F, et al. *STAR: ultrafast
   universal RNA-seq aligner.* Bioinformatics. 2013;29(1):15-21. PMID: 23104886.
   DOI: [10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)
   — used for splice-aware alignment and chimeric (fusion) read detection.

2. **fastp — read QC/trimming.** Chen S, Zhou Y, Chen Y, Gu J. *fastp: an ultra-fast
   all-in-one FASTQ preprocessor.* Bioinformatics. 2018;34(17):i884-i890. PMID: 30423086.
   DOI: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
   — adapter/quality/polyG trimming.

3. **Arriba — fusion detection.** Uhrig S, Ellermann J, Walther T, et al. *Accurate and
   efficient detection of gene fusions from RNA sequencing data.* Genome Res.
   2021;31(3):448-460. PMID: 33441414.
   DOI: [10.1101/gr.257246.119](https://doi.org/10.1101/gr.257246.119)
   — primary fusion caller; provides breakpoints and supporting read identifiers.

4. **featureCounts — expression counting.** Liao Y, Smyth GK, Shi W. *featureCounts: an
   efficient general purpose program for assigning sequence reads to genomic features.*
   Bioinformatics. 2014;30(7):923-930. PMID: 24227677.
   DOI: [10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)
   — gene-level read counting for the expression matrix.

5. **MultiQC — QC aggregation.** Ewels P, Magnusson M, Lundin S, Käller M. *MultiQC:
   summarize analysis results for multiple tools and samples in a single report.*
   Bioinformatics. 2016;32(19):3047-3048. PMID: 27312411.
   DOI: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)
   — single aggregated QC report.

6. **mosdepth — coverage.** Pedersen BS, Quinlan AR. *Mosdepth: quick coverage calculation
   for genomes and exomes.* Bioinformatics. 2018;34(5):867-868. PMID: 29096012.
   DOI: [10.1093/bioinformatics/btx699](https://doi.org/10.1093/bioinformatics/btx699)
   — target coverage and % bases at 50x/100x/200x/500x/1000x.

7. **RSeQC — RNA-seq QC.** Wang L, Wang S, Li W. *RSeQC: quality control of RNA-seq
   experiments.* Bioinformatics. 2012;28(16):2184-2185. PMID: 22743226.
   DOI: [10.1093/bioinformatics/bts356](https://doi.org/10.1093/bioinformatics/bts356)
   — read distribution and gene-body coverage (optional).

8. **SAMtools — BAM handling.** Danecek P, Bonfield JK, Liddle J, et al. *Twelve years of
   SAMtools and BCFtools.* GigaScience. 2021;10(2):giab008. PMID: 33590861.
   DOI: [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
   — BAM indexing and flagstat.

9. **Exome-capture RNA-seq for FFPE (assay rationale).** Cieslik M, Chugh R, Wu YM, et al.
   *The use of exome capture RNA-seq for highly degraded RNA with application to clinical
   cancer sequencing.* Genome Res. 2015;25(9):1372-1381. PMID: 26253700.
   DOI: [10.1101/gr.189621.115](https://doi.org/10.1101/gr.189621.115)
   — justifies capture-based RNA-seq for degraded FFPE clinical samples.

## Additional tools and standards (primary sources)

10. **STAR-Fusion — second fusion caller.** Haas BJ, Dobin A, Li B, et al. *Accuracy
    assessment of fusion transcript detection via read-mapping and de novo fusion
    transcript assembly-based methods.* Genome Biol. 2019;20:213.
    DOI: [10.1186/s13059-019-1842-9](https://doi.org/10.1186/s13059-019-1842-9)
    — orthogonal fusion caller; provides breakpoints and supporting read names.

11. **Picard toolkit.** Broad Institute. *Picard Toolkit.* GitHub:
    https://broadinstitute.github.io/picard/ — RnaSeqMetrics, HsMetrics (on-target),
    InsertSizeMetrics, BedToIntervalList, CreateSequenceDictionary.

12. **Nextflow — workflow engine.** Di Tommaso P, Chatzou M, Floden EW, et al. *Nextflow
    enables reproducible computational workflows.* Nat Biotechnol. 2017;35(4):316-319.
    DOI: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

13. **nf-core — pipeline standards.** Ewels PA, Peltzer A, Fillinger S, et al. *The nf-core
    framework for community-curated bioinformatics pipelines.* Nat Biotechnol.
    2020;38(3):276-278. DOI: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)

14. **AMP/ASCO/CAP — variant interpretation standard.** Li MM, Datto M, Duncavage EJ, et al.
    *Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in
    Cancer.* J Mol Diagn. 2017;19(1):4-23.
    DOI: [10.1016/j.jmoldx.2016.10.002](https://doi.org/10.1016/j.jmoldx.2016.10.002)
    — four-tier framework referenced for clinical reporting.
