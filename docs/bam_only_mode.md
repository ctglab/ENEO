# Running ENEO from BAM files

If your samples have already been aligned, you can skip trimming, rRNA removal, STAR
alignment, Salmon quantification, HLA typing, and neoantigen prediction, and run only
the variant calling and annotation branch. This requires three manual steps: editing
the Snakefile, preparing the sample sheet, and placing your BAM files where the
pipeline expects them.

---

## Step 1 — Edit `workflow/Snakefile`

### 1a. Remove upstream and downstream includes

Comment out every `include:` that you do not need. Keep only the rules that operate
on BAMs and VCFs:

```python
# include: "rules/index.smk"          # genome index — not needed
# include: "rules/reads_trimming.smk" # fastp + sortmerna — not needed
# include: "rules/alignment.smk"      # STAR — not needed
# include: "rules/quantification.smk" # Salmon — not needed
# include: "rules/HLA_typing.smk"     # T1K — not needed
# include: "rules/pMHC.smk"           # netMHCpan — not needed
# include: "rules/reporting.smk"      # MultiQC — not needed

include: "rules/bam_cleaning.smk"
include: "rules/base_recalibration.smk"
include: "rules/strelka.smk"
include: "rules/deepvariant.smk"
include: "rules/filter_calls.smk"
include: "rules/annotate_variants.smk"
```

### 1b. Simplify `rule targets`

Replace the existing `rule targets` with one that only requests the final annotated
VCF per patient:

```python
rule targets:
    input:
        expand(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["VCF_out"],
                "{patient}_final_passonly.vcf.gz.tbi"
            ),
            patient=patients,
        ),
```

---

## Step 2 — Prepare `units.csv`

The pipeline reads sample metadata from `units.csv`. In normal mode this file has
columns `patient`, `fq1`, and `fq2`. In BAM-only mode those columns are never
accessed (the rules that call `get_fastq()` are no longer included), but the file
must still be parseable with a `patient` column so that the wildcard constraints in
`common.smk` resolve correctly.

The safest approach is to keep the header and fill `fq1`/`fq2` with placeholder
values:

```csv
patient,fq1,fq2
pat25,NA,NA
```

The placeholder strings are never read at runtime.

Also make sure `patients.csv` lists every patient you want to process:

```csv
patient
pat25
```

---

## Step 3 — Place your BAM at the expected input path

The first BAM-processing rule (`AddGrp` in `bam_cleaning.smk`) looks for its input
at:

```
{OUTPUT_FOLDER}/mapped_reads/{patient}_Aligned.sortedByCoord.out.bam
```

where `OUTPUT_FOLDER` is the value set in `config/config_main.yaml`.

Create the `mapped_reads` directory and symlink (or copy) your BAM there with the
correct name. For a patient called `pat25`:

```bash
mkdir -p <OUTPUT_FOLDER>/mapped_reads

# symlink (preferred — no extra disk usage)
ln -s /absolute/path/to/pat25.bam \
      <OUTPUT_FOLDER>/mapped_reads/pat25_Aligned.sortedByCoord.out.bam
```

**Requirements for the input BAM:**

- Must be **coordinate-sorted**. The pipeline runs `samtools sort` internally after
  `AddGrp`, so a coordinate-sorted input is required.
- Does **not** need a pre-existing `.bai` index. `samtools_index` is run later in the
  pipeline.
- Read groups will be overwritten by `AddGrp` regardless of what is already present in
  the BAM header.

---

## Step 4 — Run the workflow

No other configuration changes are needed. Launch Snakemake as usual:

```bash
snakemake --profile workflow/profile/slurm_profile
```

The DAG will start at `AddGrp` and produce the final per-patient VCF:

```
{OUTPUT_FOLDER}/VCF_out/{patient}_final_passonly.vcf.gz
{OUTPUT_FOLDER}/VCF_out/{patient}_final_passonly.vcf.gz.tbi
```

---

## Restoring the full workflow

To go back to a full run, simply uncomment the `include:` lines and restore the
original `rule targets`. No rule files were modified, so no further changes are
needed.
