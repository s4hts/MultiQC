---
title: PURPLE
description: >
  A purity, ploidy and copy number estimator for whole genome tumor data
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/purple/purple.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

:::note
A purity, ploidy and copy number estimator for whole genome tumor data

[https://github.com/hartwigmedical/hmftools/](https://github.com/hartwigmedical/hmftools/)
:::

PURPLE combines B-allele frequency (BAF), read depth ratios, somatic variants and
structural variant breakpoints to estimate the purity and copy number profile
of a tumor sample, and also predicts gender, the MSI status, tumor mutational
load and burden, clonality and the whole genome duplication status.

### File search patterns

```yaml
purple/purity:
  fn: "*.purple.purity.tsv"
purple/qc:
  fn: "*.purple.qc"
```