# R言語による生物情報学ドキュメント例

#ドキュメント #R #RMarkdown #roxygen2 #実践例 #生物情報学 #RNA-seq #FAIR

> [!NOTE]
> このファイルは[[Composer_ProjectRulesTDD]]から分割された、ドキュメント作成の具体例を示すものです。
> 関連ガイド: [[Composer_ProjectRulesTDD]] | [[Composer_ProjectRules_Cursor実装例]]

## 1. roxygen2ドキュメント例

関数やクラスのドキュメントは以下のように記述します。これにより自動的にRパッケージのドキュメントが生成されます。

```r
#' @title Calculate the mean of a numeric vector
#' @description
#' Calculates the arithmetic mean of values in a numeric vector.
#' Missing values can be removed before calculation.
#' 
#' @param x A numeric vector whose mean is to be calculated
#' @param na.rm Logical. Should missing values be removed? Default is TRUE
#' 
#' @return The arithmetic mean of the values in \code{x}
#' 
#' @examples
#' calculate_mean(c(1, 2, 3, 4, 5))
#' calculate_mean(c(1, 2, NA, 4, 5), na.rm = TRUE)
#' 
#' @export
calculate_mean <- function(x, na.rm = TRUE) {
  if (!is.numeric(x)) {
    stop("Input must be numeric")
  }
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  mean(x, na.rm = na.rm)
}
```

## 2. R Markdownドキュメント例

RNA-seq解析のステップをR Markdownで記録する例です。これによりコードと解説が一体化した再現可能なドキュメントが作成できます。

```r
---
title: "DESeq2によるRNA-seq差分発現解析"
author: "解析者名"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: show
    fig_caption: true
    self_contained: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 10,
  fig.height = 7,
  fig.path = "figures/",
  dev = c("png", "pdf")
)

# 環境管理
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

# ライブラリ読み込み
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(futile.logger)
library(cli)

# ロギング設定
log_file <- paste0("analysis_", format(Sys.time(), "%Y%m%d_%H%M"), ".log")
flog.appender(appender.file(log_file))
flog.threshold(INFO)
```

## 1. 研究背景と目的

[解析の背景と目的を説明。具体的な仮説や検証したいことを明記]

## 2. データ概要

```{r data_overview}
# メタデータ読み込み
sample_info <- read_csv("metadata/sample_info.csv")
knitr::kable(sample_info)

# データ読み込み
counts <- read.csv("data/raw_counts.csv", row.names = 1)
```

## 3. 差分発現解析

[解析コード、結果、可視化]

## 4. Gene Ontology解析

[GO解析の詳細]

## 5. 結論と次のステップ

[結論と今後の方向性]

```{r session_info}
sessionInfo()
```
```

## 3. ISA-Tab形式メタデータ例

ISA-Tabは実験メタデータを標準化された形式で記述するためのフォーマットです。以下は簡略化した例です。

### investigation.tsv (研究全体のメタデータ)

```
ONTOLOGY SOURCE REFERENCE
Term Source Name	OBI	EFO	NCBI Taxonomy
Term Source File	http://purl.obolibrary.org/obo/obi.owl	http://www.ebi.ac.uk/efo/	http://www.ncbi.nlm.nih.gov/taxonomy

INVESTIGATION
Investigation Identifier	RNA-Seq-Exp-2023
Investigation Title	Gene expression profiling of treatment X on cell line Y
Investigation Description	This study examines differential gene expression in response to treatment X on cell line Y
Investigation Submission Date	2023-10-15
Investigation Public Release Date	2024-04-15

INVESTIGATION PUBLICATIONS
Investigation PubMed ID	
Investigation Publication DOI	
Investigation Publication Author List	
Investigation Publication Title	
Investigation Publication Status	Planned

INVESTIGATION CONTACTS
Investigation Person Last Name	Smith
Investigation Person First Name	John
Investigation Person Email	john.smith@example.com
Investigation Person Affiliation	University of Science
Investigation Person Roles	corresponding author
```

### study.tsv (実験デザインのメタデータ)

```
STUDY
Study Identifier	RNA-Seq-Study-01
Study Title	Treatment X effects on gene expression
Study Description	Examination of transcriptional changes after treatment X exposure
Study Submission Date	2023-10-15
Study Public Release Date	2024-04-15
Study File Name	s_rna_seq_study.tsv

STUDY DESIGN DESCRIPTORS
Study Design Type	case-control design
Study Design Type Term Accession Number	OBI:0500013
Study Design Type Term Source REF	OBI

STUDY FACTORS
Study Factor Name	treatment	exposure_time
Study Factor Type	chemical compound	time

STUDY ASSAYS
Study Assay Measurement Type	transcription profiling
Study Assay Measurement Type Term Accession Number	OBI:0000424
Study Assay Measurement Type Term Source REF	OBI
Study Assay Technology Type	RNA sequencing
Study Assay Technology Type Term Accession Number	OBI:0001271
Study Assay Technology Type Term Source REF	OBI
Study Assay Technology Platform	Illumina NovaSeq 6000
Study Assay File Name	a_rna_seq_assay.tsv
```

### assay.tsv (実験手法と測定のメタデータ)

```
Sample Name	Protocol REF	Extract Name	Protocol REF	Library Name	Protocol REF	Sequencing Assay Name	Raw Data File
sample1	RNA extraction	extract1	RNA-seq library prep	lib1	RNA sequencing	seq1	sample1_R1.fastq.gz;sample1_R2.fastq.gz
sample2	RNA extraction	extract2	RNA-seq library prep	lib2	RNA sequencing	seq2	sample2_R1.fastq.gz;sample2_R2.fastq.gz
sample3	RNA extraction	extract3	RNA-seq library prep	lib3	RNA sequencing	seq3	sample3_R1.fastq.gz;sample3_R2.fastq.gz
sample4	RNA extraction	extract4	RNA-seq library prep	lib4	RNA sequencing	seq4	sample4_R1.fastq.gz;sample4_R2.fastq.gz
```

## 4. 環境情報の記録例

解析環境の記録関数の例です。これは再現性確保のために重要です。

```r
# 環境情報記録の関数
log_session_info <- function(file_prefix = "session_info") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file_name <- paste0(file_prefix, "_", timestamp, ".txt")
  
  sink(file_name)
  cat("## R Session Info\n\n")
  print(sessionInfo())
  cat("\n\n## Installed Packages\n\n")
  print(installed.packages()[, c("Package", "Version")])
  cat("\n\n## System Information\n\n")
  cat("Date: ", as.character(Sys.Date()), "\n")
  cat("System: ", Sys.info()["sysname"], "\n")
  cat("User: ", Sys.info()["user"], "\n")
  sink()
  
  return(file_name)
}

# 解析開始時に呼び出す
session_file <- log_session_info("rnaseq_differential_expression")
```

## 5. BioComputeObject (BCO)例

以下はJSON形式のBCOの簡略化した例です。実際のBCOはより詳細なメタデータを含みます。

```json
{
  "object_id": "https://example.org/compute_objects/RNA-seq-analysis-1234",
  "spec_version": "https://w3id.org/ieee/ieee-2791-schema/2791object.json",
  "etag": "abcdef123456",
  "provenance_domain": {
    "name": "RNA-seq differential expression workflow",
    "version": "1.0.0",
    "created": "2023-10-15T12:00:00Z",
    "modified": "2023-10-15T12:00:00Z",
    "contributors": [
      {
        "name": "John Smith",
        "email": "john.smith@example.com",
        "affiliation": "University of Science",
        "contribution": ["authoredBy"]
      }
    ],
    "license": "https://spdx.org/licenses/CC-BY-4.0.html"
  },
  "usability_domain": [
    "This workflow performs differential expression analysis on RNA-seq data",
    "Designed for paired-end Illumina RNA-seq data from human samples",
    "Performs read QC, alignment, quantification and differential expression with DESeq2"
  ],
  "description_domain": {
    "keywords": ["RNA-seq", "differential expression", "transcriptomics"],
    "platform": ["Illumina NovaSeq"],
    "pipeline_steps": [
      {
        "step_number": 1,
        "name": "Quality Control",
        "description": "FastQC and MultiQC for read quality assessment",
        "input_list": [
          {
            "uri": "file://data/raw/*.fastq.gz",
            "access_time": "2023-10-01T10:00:00Z"
          }
        ],
        "output_list": [
          {
            "uri": "file://results/qc/multiqc_report.html"
          }
        ]
      },
      {
        "step_number": 2,
        "name": "Read Alignment",
        "description": "STAR aligner for mapping reads to reference genome",
        "prerequisites": ["Quality Control"],
        "input_list": [
          {
            "uri": "file://data/raw/*.fastq.gz"
          },
          {
            "uri": "file://reference/genome.fa"
          }
        ],
        "output_list": [
          {
            "uri": "file://results/aligned/*.bam"
          }
        ]
      }
    ]
  },
  "execution_domain": {
    "script": [
      {
        "uri": "file://code/01_quality_control.Rmd"
      },
      {
        "uri": "file://code/02_alignment.sh"
      },
      {
        "uri": "file://code/03_quantification.R"
      },
      {
        "uri": "file://code/04_differential_expression.Rmd"
      }
    ],
    "script_driver": "Rscript",
    "software_prerequisites": [
      {
        "name": "R",
        "version": "4.2.2",
        "uri": "https://www.r-project.org/"
      },
      {
        "name": "DESeq2",
        "version": "1.38.0",
        "uri": "https://bioconductor.org/packages/DESeq2/"
      }
    ],
    "external_data_endpoints": [],
    "environment_variables": {}
  },
  "parametric_domain": [
    {
      "param": "p_value_threshold",
      "value": "0.05",
      "step": "Differential Expression"
    },
    {
      "param": "fold_change_threshold",
      "value": "1.5",
      "step": "Differential Expression"
    }
  ],
  "io_domain": {
    "input_subdomain": [
      {
        "name": "Raw RNA-seq Reads",
        "uri": "file://data/raw/*.fastq.gz",
        "access_time": "2023-10-01T10:00:00Z"
      },
      {
        "name": "Sample Metadata",
        "uri": "file://data/metadata/sample_info.csv",
        "access_time": "2023-10-01T10:00:00Z"
      }
    ],
    "output_subdomain": [
      {
        "name": "Differential Expression Results",
        "uri": "file://results/tables/differential_expression.csv",
        "access_time": "2023-10-15T12:00:00Z"
      },
      {
        "name": "Visualization",
        "uri": "file://results/figures/volcano_plot.pdf",
        "access_time": "2023-10-15T12:00:00Z"
      }
    ]
  }
}
```

## 6. Snakemakeワークフロー例

パイプラインをSnakemakeで記述した例です。これにより解析の再現性と自動化が可能になります。

```python
# Snakefile for RNA-seq analysis

# Configuration
configfile: "workflow/config.yaml"

# Global variables
SAMPLES = config["samples"]
REFERENCE_GENOME = config["reference_genome"]
OUTPUT_DIR = config["output_dir"]

# All target rule
rule all:
    input:
        # Quality control
        "results/qc/multiqc_report.html",
        # Alignment results
        expand("results/aligned/{sample}.bam", sample=SAMPLES),
        # Quantification results
        "results/counts/raw_counts.csv",
        # Differential expression
        "results/diffexp/deseq2_results.csv",
        # Final report
        "results/reports/complete_analysis.html"

# Quality Control
rule fastqc:
    input:
        "data/raw/{sample}_R{read}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}_R{read}_fastqc.html",
        zip="results/qc/fastqc/{sample}_R{read}_fastqc.zip"
    log:
        "logs/fastqc/{sample}_R{read}.log"
    threads: 2
    shell:
        "fastqc -o results/qc/fastqc -t {threads} {input} > {log} 2>&1"

rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_R{read}_fastqc.zip", 
               sample=SAMPLES, read=[1, 2])
    output:
        "results/qc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        "multiqc results/qc/fastqc -o results/qc > {log} 2>&1"

# Alignment
rule star_align:
    input:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
        ref=REFERENCE_GENOME
    output:
        bam="results/aligned/{sample}.bam"
    log:
        "logs/star/{sample}.log"
    threads: 12
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.ref} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand gunzip -c \
             --outSAMtype BAM SortedByCoordinate \
             --outBAMsortingThreadN 6 \
             --outFileNamePrefix results/aligned/{wildcards.sample}_ \
             --quantMode GeneCounts > {log} 2>&1
        
        mv results/aligned/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        """

# Quantification
rule featureCounts:
    input:
        bams=expand("results/aligned/{sample}.bam", sample=SAMPLES),
        annotation=config["annotation"]
    output:
        counts="results/counts/raw_counts.csv"
    log:
        "logs/featurecounts/featurecounts.log"
    threads: 8
    shell:
        """
        featureCounts -T {threads} \
                     -p \
                     -t exon \
                     -g gene_id \
                     -a {input.annotation} \
                     -o results/counts/counts.txt \
                     {input.bams} > {log} 2>&1
        
        # Convert to CSV format
        Rscript -e "
        counts <- read.table('results/counts/counts.txt', header=TRUE, skip=1)
        counts <- counts[,c(1, 7:ncol(counts))]
        colnames(counts)[1] <- 'gene_id'
        colnames(counts)[2:ncol(counts)] <- sub('.*aligned\\\/', '', sub('.bam', '', colnames(counts)[2:ncol(counts)]))
        write.csv(counts, '{output.counts}', row.names=FALSE)
        " >> {log} 2>&1
        """

# Differential Expression
rule deseq2:
    input:
        counts="results/counts/raw_counts.csv",
        metadata=config["metadata"]
    output:
        results="results/diffexp/deseq2_results.csv",
        plots="results/diffexp/figures/volcano_plot.pdf"
    log:
        "logs/deseq2/deseq2.log"
    script:
        "scripts/run_deseq2.R"

# Generate Final Report
rule final_report:
    input:
        qc="results/qc/multiqc_report.html",
        counts="results/counts/raw_counts.csv",
        diffexp="results/diffexp/deseq2_results.csv"
    output:
        report="results/reports/complete_analysis.html"
    log:
        "logs/report/final_report.log"
    shell:
        """
        Rscript -e "
        rmarkdown::render('code/complete_analysis.Rmd', 
                          output_file='{output.report}',
                          params=list(
                            qc_report='{input.qc}',
                            counts='{input.counts}',
                            diffexp='{input.diffexp}'
                          ))
        " > {log} 2>&1
        """
```

## 7. ディレクトリ構造例

RNA-seqプロジェクトの標準的なディレクトリ構造です。この構造に従うことで、データの整理と管理が容易になります。

```
project_root/
├── README.md                  # プロジェクト概要
├── .gitignore                 # Git除外ファイル設定
├── .dvcignore                 # DVC除外ファイル設定
├── renv/                      # R環境管理
├── data/                      # データファイル（多くはgitignore対象）
│   ├── raw/                   # 生データファイル
│   ├── processed/             # 処理済みデータ
│   └── metadata/              # メタデータファイル
├── code/                      # 解析スクリプト
│   ├── 01_quality_control.Rmd # 品質管理
│   ├── 02_normalization.Rmd   # 正規化
│   ├── 03_differential_exp.Rmd # 差分発現解析
│   └── utils/                 # 共通ユーティリティ関数
├── workflow/                  # ワークフロー定義ファイル
│   ├── Snakefile              # Snakemakeワークフロー
│   └── config.yaml            # 設定ファイル
├── results/                   # 解析結果
│   ├── figures/               # 生成された図表
│   ├── tables/                # 生成された表
│   └── reports/               # レポートHTML/PDF
└── docs/                      # プロジェクトドキュメント
    ├── biocompute_objects/    # BCOファイル
    └── protocol/              # 実験・解析プロトコル
``` 