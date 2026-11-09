# Software versions (RNA-seq)

| Tool | Version | Used in | Source of version |
|------|---------|---------|-------------------|
| Cutadapt | 5.0 | step 00 trimming | MultiQC report |
| FastQC | 0.12.1 | step 00 QC | MultiQC report |
| UMI-tools | 1.1.6 | step 00 (SMART-seq UMI extraction/deduplication) | MultiQC report |
| STAR | 2.7.4a (rDNA index) | steps 00, 01, 02b–02d | genomeParameters.txt |
| Salmon | 1.10.3 | step 00 (SMART-seq), 03c | MultiQC report |
| featureCounts (Subread) | 2.1.1 | steps 03a, 03b | featureCounts output header |
| bowtie2 | TBD | step 02a | |
| Samtools | TBD | steps 00, 02 | |
| deepTools | TBD | steps 00, 02e | |
