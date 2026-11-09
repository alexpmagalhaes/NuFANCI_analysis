# NuFANCI: analysis code

> Part of the manuscript **"Dual and latent condensate localization signatures encoded in intrinsically disordered protein regions"**.
> Author and contact: **Alexandre P. Magalhaes**, Email: magalhae@molgen.mpg.de

Analysis scripts for the NuFANCI study: proteomics, DNA-seq, RNA-seq and nanopore (ONT) sequencing of nucleoli, nuclear speckles and related experiments (for example Actinomycin D treatment).

Raw and processed sequencing data (FASTQ, BAM, BED, bigWig, count matrices) and MS raw files are **not** stored here. They are deposited in public repositories:

| Data | Repository | Accession |
|------|------------|-----------|
| DNA-seq | GEO | [GSE353840](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353840) |
| RNA-seq | GEO | [GSE335756](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335756) |
| ONT (nanopore, CpG methylation) | GEO | GSE_TBD |
| MS proteomics | PRIDE | [PXD069855](https://www.ebi.ac.uk/pride/archive/projects/PXD069855) |

## Related software

- **[FPCaller](https://github.com/alexpmagalhaes/FPCaller)**: the domain caller used for the DNA-seq analysis (NAD/SAD calling). It is released as a separate Python package.

## Repository layout

```
./
├── NuFANCI_DNAseq/          # DNA-seq
├── NuFANCI_RNAseq/          # RNA-seq
├── NuFANCI_MS_Proteomics/   # mass-spectrometry proteomics
└── NuFANCI_ONT/             # Oxford Nanopore sequencing
```

Each module has the same structure:

```
<module>/
├── README.md      # what the module does and the order to run the scripts
├── scripts/       # numbered in run order: 01_*, 02_*, ...
├── results/       # small result tables and figures (no large files)
└── envs/          # software versions / conda environment files
```

## How to reproduce

1. Download the data from the accessions above.
2. Set up the environment from `<module>/envs/`.
3. Run the scripts in `<module>/scripts/` in numeric order (see each module's README).

## Citation

This code is part of the manuscript:

> Dual and latent condensate localization signatures encoded in intrinsically disordered protein regions.
> (Journal, year and DOI: TBD)

## Contact

**Alexandre P. Magalhaes**, Email: magalhae@molgen.mpg.de
