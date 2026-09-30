# Bioinformatics Teaching Lab

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/YOUR-USERNAME/bioinfo-teaching-lab?quickstart=1)

A ready to use cloud workspace for learning the modern bioinformatics stack. Everything runs in your web browser through GitHub Codespaces, so there is nothing to install on your own computer.

## What is inside

| Area | Tools |
| --- | --- |
| Workflow managers | Nextflow, nf-core tools, nf-test, Snakemake |
| Software environments | conda and mamba (Miniforge), Docker |
| Read QC and trimming | FastQC, fastp, MultiQC, SeqKit |
| Alignment and quantification | BWA, minimap2, Salmon |
| BAM, VCF and intervals | samtools, bcftools, htslib (bgzip, tabix), bedtools |
| Python analysis | pandas, Biopython, matplotlib, seaborn, Jupyter notebooks |

You also get a tiny simulated dataset (3 samples, a 50 kb genome, 40 planted SNPs) and a small Nextflow pipeline that goes from raw reads all the way to variant calls.

## Getting started

1. Near the top of this page, click **Use this template**, then **Create a new repository**. This gives you your own copy where your work is saved.
2. In your new repository, click **Code**, open the **Codespaces** tab, and click **Create codespace on main**.
3. Wait for setup to finish. The very first start can take several minutes. Restarting it later is much faster.
4. When the terminal prints **Setup complete**, run:

```bash
nextflow run main.nf
```

If everything worked you now have a `results/` folder. Keep reading to understand what just happened.

## Look after your free hours

Every personal GitHub account gets a free monthly allowance of Codespaces time and storage. On the default 2 core machine that is roughly 60 hours per month. Verified students with the GitHub Student Developer Pack get more.

* Always choose the **2 core** machine. Bigger machines use up your hours faster.
* **Stop** your codespace when you finish. Go to [github.com/codespaces](https://github.com/codespaces), click the three dots next to it, and choose **Stop codespace**. Idle codespaces stop on their own after 30 minutes, but those minutes still count.
* **Delete** codespaces you no longer need. A stopped codespace still uses storage.
* Commit and push your work before deleting, otherwise it is gone.

## What is in this repository

```
.
├── .devcontainer/
│   ├── devcontainer.json       How the codespace is set up (tools, extensions, machine size)
│   ├── Dockerfile              Builds the image: Ubuntu + conda + the bioinfo environment
│   └── post-create.sh          Runs once at startup: makes the mock data, checks tools
├── .github/workflows/
│   └── build-image.yml         Optional: prebuilds the image for faster startups
├── conf/
│   └── codespaces.config       Resource limits for running nf-core pipelines here
├── data/
│   └── generate_mock_data.py   Simulates the reads, reference and truth set
├── scripts/
│   └── check_env.sh            Checks that every tool is installed
├── tests/
│   └── main.nf.test            An nf-test test for the pipeline
├── environment.yml             The conda environment (the full tool list)
├── main.nf                     The Nextflow pipeline
├── nextflow.config             Pipeline settings and profiles
└── nf-test.config              nf-test settings
```

The mock data lives in `data/mock/`. It is generated fresh in every codespace, so it is not stored in git.

## Lessons

### Lesson 1: Meet the data

```bash
ls -lh data/mock
seqkit stats data/mock/*.fastq.gz
zcat data/mock/sampleA_R1.fastq.gz | head -8
grep ">" data/mock/reference.fa
column -t data/mock/truth_snps.tsv | head
```

Questions to answer:

* How many read pairs does each sample have, and how long are the reads?
* What does each of the 4 lines of a FASTQ record mean?
* `truth_snps.tsv` lists the SNPs hidden in the data. Later you will check whether the pipeline finds them.

### Lesson 2: Quality control by hand

```bash
mkdir -p scratch && cd scratch
fastqc ../data/mock/sampleA_R1.fastq.gz ../data/mock/sampleC_R1.fastq.gz -o .
fastp -i ../data/mock/sampleC_R1.fastq.gz -I ../data/mock/sampleC_R2.fastq.gz \
      -o sampleC_R1.trim.fastq.gz -O sampleC_R2.trim.fastq.gz \
      --detect_adapter_for_pe --html sampleC_fastp.html
cd ..
```

Open the HTML reports (see "Viewing HTML reports" below) and compare sampleA with sampleC. What is wrong with sampleC, and what did fastp do about it?

### Lesson 3: Run the pipeline

```bash
nextflow run main.nf
```

Open `main.nf` next to the terminal and find each step:

* **Processes** (FASTQC, FASTP, BWA_MEM and so on) each wrap one tool.
* **Channels** carry data from one process to the next.
* The **workflow** block at the bottom wires everything together.

Nextflow runs every task in its own folder under `work/`. The final files are copied to `results/`. The full log is in `.nextflow.log`.

### Lesson 4: Explore the results

* `results/multiqc/multiqc_report.html`: one report combining FastQC, fastp, samtools and bcftools
* `results/pipeline_info/report.html`: how long each task took and how much memory it used
* `results/variants/cohort.vcf.gz`: the variant calls

List the calls in a readable table:

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' results/variants/cohort.vcf.gz
```

Compare them with `data/mock/truth_snps.tsv`. Did the pipeline find every planted SNP? Are the genotypes right? Did it report anything that is not in the truth set?

### Lesson 5: Resume, parameters and changes

```bash
nextflow run main.nf -resume                        # nothing reruns, everything is cached
nextflow run main.nf -resume --outdir results_v2    # change a parameter
nextflow log                                        # history of your runs
```

Now edit a process in `main.nf` (for example add `--cut_tail` to the fastp command) and run with `-resume` again. Which steps rerun, and why?

### Lesson 6: Same pipeline, different software sources

```bash
nextflow run main.nf -profile conda  --outdir results_conda
nextflow run main.nf -profile docker --outdir results_docker
```

* With `-profile conda`, Nextflow builds a small conda environment for each step, using the exact versions listed in `nextflow.config`.
* With `-profile docker`, each step runs inside its own container image from BioContainers.

Discuss: why do real pipelines pin exact tool versions, and why are containers the most common choice for sharing pipelines?

### Lesson 7: Test your pipeline

```bash
nf-test test tests/main.nf.test
```

Open `tests/main.nf.test`. Try adding a new assertion, for example that `results/alignment/sampleA.sorted.bam` exists.

### Lesson 8: Run a real nf-core pipeline

```bash
nf-core pipelines list
nextflow run nf-core/demo -profile test,docker -c conf/codespaces.config --outdir results_nfcore
```

`conf/codespaces.config` caps every task at 2 CPUs and 6 GB of memory so community pipelines fit on this machine.

### Make new data

Change the depth, the random seed, or the number of SNPs and see how the results change:

```bash
python data/generate_mock_data.py --outdir data/mock --coverage 10 --seed 7
nextflow run main.nf --outdir results_lowcov
```

## Viewing HTML reports

Pick whichever you prefer:

* In the file explorer, right click an `.html` file and choose **Show Preview**.
* Or serve the results folder and open the link Codespaces gives you:

```bash
python -m http.server 8000 --directory results
```

* Or right click the file and choose **Download**.

## Adding more tools

Install into the running environment for a quick try:

```bash
mamba install -n bioinfo <tool-name>
```

To make it permanent for everyone, add the tool to `environment.yml`, commit, and rebuild: open the Command Palette (Ctrl+Shift+P or Cmd+Shift+P) and run **Codespaces: Rebuild Container**.

## Troubleshooting

| Problem | Fix |
| --- | --- |
| `command not found` | Run `conda activate bioinfo`, or open a new terminal |
| A pipeline step failed | Read the error, then look at `.command.err` inside the task folder under `work/` that Nextflow printed |
| "Process requirement exceeds available CPUs" | Add `-c conf/codespaces.config` to your command |
| Disk is full | `rm -rf work results_*` and `docker system prune -af` |
| Everything is broken | Delete the codespace and create a new one. Your pushed code is safe on GitHub |

Check that all tools are present at any time:

```bash
bash scripts/check_env.sh
```

## For instructors

**Faster startups.** By default each codespace builds the image from the Dockerfile, which takes several minutes. To make startups much faster, prebuild the image once:

1. Go to the **Actions** tab, choose **Build codespace image**, and click **Run workflow**.
2. When it finishes, open your profile's **Packages** page, open the new package, go to **Package settings**, and set visibility to **Public**.
3. In `.devcontainer/devcontainer.json`, replace the whole `"build": { ... }` block with:

```json
"image": "ghcr.io/YOUR-USERNAME/bioinfo-teaching-lab:latest",
```

After that, the workflow rebuilds the image automatically whenever you change `environment.yml` or the Dockerfile.

**Resetting a run.** `rm -rf work results* .nextflow* .nf-test` brings the workspace back to a clean state.
