# In silico detection of *Xylella fastidiosa* transcripts in publicly available insect vector RNA-seq datasets

## Aim

The aim of this project is to develop a reproducible bioinformatic
pipeline in Python to mine publicly available insect vector RNA-seq
datasets for evidence of *Xylella fastidiosa* infection and recover
bacterial transcripts captured incidentally during host-focused
sequencing experiments.

------------------------------------------------------------------------

## Background

*Xylella fastidiosa* is a xylem-limited bacterial plant pathogen
transmitted by xylophagous insect vectors such as sharpshooters and
spittlebugs. While its transcriptional landscape within plant hosts has
been studied to some extent, the bacterial gene expression profile
during colonization of the insect foregut remains largely unexplored.
This represents a critical knowledge gap, as successful transmission
depends on bacterial adhesion, biofilm formation, and persistence within
the vector's feeding apparatus.

Currently, there are very few (if any) RNA-seq datasets specifically
designed to capture the transcriptome of *X. fastidiosa* within insect
vectors. However, multiple transcriptomic studies focusing on insect
physiology, immunity, or feeding behavior may have inadvertently
sequenced bacterial RNA present during natural or experimental
infections.

The use of such datasets presents an opportunity to computationally
detect and reconstruct pathogen-derived transcripts through post hoc
analysis of raw sequencing reads. This approach resembles dual RNA-seq
or metatranscriptomic profiling strategies and can provide preliminary
insight into bacterial gene expression during vector colonization
without the need for dedicated infection experiments.

------------------------------------------------------------------------

## Project Objectives

This project aims to:

-   Determine whether publicly available insect RNA-seq datasets contain
    detectable levels of *X. fastidiosa* transcripts.
-   Develop a Python-based pipeline capable of identifying and
    extracting pathogen-derived reads from host-centric sequencing
    datasets.
-   Reconstruct bacterial transcripts and perform preliminary functional
    annotation of expressed genes.
-   Explore the feasibility of using host-focused RNA-seq datasets for
    pathogen transcriptomic profiling.

------------------------------------------------------------------------

## Bioinformatic Pipeline Overview

The proposed pipeline will consist of the following steps:

1.  Raw RNA-seq reads will be retrieved from public repositories such as
    the NCBI Sequence Read Archive (SRA).
2.  Quality control and adapter trimming will be performed using
    standard preprocessing tools.
3.  Reads will be mapped against a reference genome of *Xylella
    fastidiosa* to identify pathogen-derived sequences.
4.  Mapped reads will be extracted and assembled into putative bacterial
    transcripts.
5.  Recovered transcripts will be annotated using available gene models
    to identify expressed virulence-associated functions such as
    adhesins or biofilm-related proteins.

Custom Python scripts will be developed to automate read filtering,
mapping statistics parsing, transcript extraction, and downstream data
formatting.

------------------------------------------------------------------------

## Expected Outcomes

This project is expected to result in:

-   A functional and modular Python pipeline for pathogen transcript
    detection in host RNA-seq datasets.
-   Identification of candidate *X. fastidiosa* transcripts expressed
    during vector colonization.
-   A proof-of-concept demonstration of indirect pathogen
    transcriptomics through public data mining.

------------------------------------------------------------------------

## Project update 04/13/2026

### Status Summary
The current phase of the project has focused heavily on developing a robust, HPC-ready Preprocessing Pipeline (preprocess_homalodisca.sh). While the initial scope intended to reach deeper into the downstream analysis (xylella_miner.py), the technical challenges and the necessity for a high-quality data foundation shifted our primary efforts toward the preprocessing stage.

 ### Key Milestones Achieved:
- Infrastructure & Scalability: Transitioned the host-filtering engine from BWA to STAR, optimized specifically for the Alabama HPC environment (Slurm).

- Modularity & Reproducibility: Implemented automated genome indexing, dependency validation, and flexible argument parsing. This ensures that the pipeline can be reused for any host-pathogen study beyond this specific project.

- Data Integrity: Developed rigorous checks for paired-end synchronization and host-read removal, ensuring that the input for the next stage is clean and reliable.

### Technical Debt & Next Steps
- Xylella Miner Development: Due to the extensive time allocated to ensuring the reproducibility and efficiency of the preprocessing script, the xylella_miner.py tool is currently in its early conceptual/boilerplate stage.

- Justification: We prioritized a solid preprocessing foundation. In bioinformatics, downstream results are only as good as the input data ("garbage in, garbage out"). By delivering a modular and robust preprocessing script, we have established a reproducible framework that guarantees the project's long-term viability and ease of scaling.
