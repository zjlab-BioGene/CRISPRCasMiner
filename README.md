# CRISPRCasMiner v.1.5.0
A CRISPR-Cas systems mining pipeline.[Colab_Notebook](https://colab.research.google.com/drive/1PYo_vFefUnPWgFLQ5q3Oxu2pTtx9BvzY?usp=sharing)

`Input`: Metagenome-assembled genomes/contigs.

`Output`: Known CRISPR-Cas systems & suspicious proteins adjacent to CRISPR arrays.

## Authors

| Author | Email |
| ------ | ----- |
| Wenhui Li | liwh@zhejianglab.com, liwh@tongji.edu.cn |
| Runze Cai | cairz@zhejianglab.com |
| Wuke Wang | wangwk@zhejianglab.com |

## Installation

1. download CRISPRCasMiner from github

```shell
git clone https://github.com/zjlab-BioGene/CRISPRCasMiner.git
```

2. unpack database
```shell
tar -xvzf CRISPRCasMiner/data/Profiles.tar.gz
mv Profiles/ CRISPRCasMiner/data/
rm CRISPRCasMiner/data/Profiles.tar.gz
```

3. environment

Create environment.
```shell
conda create -n ccminer python=3.10
conda activate ccminer
```

Install dependencies for cctyper.
```shell
conda install -c bioconda hmmer -y
conda install -c bioconda prodigal -y
conda install -c bioconda minced -y
conda install -c bioconda bioawk -y
python -m pip install 'drawSvg~=1.9'
python -m pip install 'cctyper~=1.8.0'
```

Install dependencies for ccminer.
```shell
conda install -c bioconda mafft -y
conda install -c bioconda viennarna -y
conda install -c conda-forge webencodings -y
conda install -c bioconda seqkit -y
conda install -c bioconda seqtk -y
```

## Run

1. Run `cctyper` and `prodigal`.
First, run `cctyper`:
```shell
cd CRISPRCasMiner && rm -rf output && mkdir output
cctyper example/input_test.fna output/01_cctyper \
    --db data \
    --prodigal meta \
    --keep_tmp
```
Then, run `prodigal`:
```shell
prodigal -i example/input_test.fna -d output/01_cctyper/genes.cds -p meta > /dev/null
```

2. Run ccminer.

See parameters of ccminer:
```shell
python ccminer/ccminer.py 
```
or 
```shell
python ccminer/ccminer.py -h
```
Then, run ccminer:
```shell
python ccminer/ccminer.py example/input_test.fna output/02_ccminer \
    --cctyper_path output/01_cctyper \
    --database_name my_project_name \
    --name my_sample_name \
    --db data \
    --prodigal meta \
    --span 10 \
    --keep_tmp
```
Fetch suspicous proteins:
```shell
awk '$12=="primary" {printf ">%s\n%s\n", $9,$13}' output/02_ccminer/out.tab > output/suspicious.faa
```
Finally, have a look at the proteins:
```shell
cat output/suspicious.faa
```
---------------------------------------------------------------

## *Note*

Two bioinformatics pipelines were executed on your metagenome-assembled data:

- `cctyper`. All known CRISPR-Cas systems with certainty.

- `ccminer`. All suspious proteins located near the CRISPR arrays, excluding the already known systems.


The folder tree of output directory:

```
YOUR_PROJECT_output
output
├── 01_cctyper              # output of cctyper v1.8.0
├── 02_ccminer              # output of ccminer v1.5.0
│   ├── arguments.tab
│   ├── crisprs_calib.tab
│   ├── genes.tab
│   ├── minced.out
│   ├── nearCRISPR.txt
│  *├── out.tab             # proteins adjacent to CRISPRs
│   └── seqLen.tab
└── suspicious.faa          # suspicious Cas proteins
```

If no known CRISPR-Cas systems were identified by cctyper, then pay your attention to the suspious proteins adjacent to the CRISPR array. For the output table of ccminer, `out.tab` presents the following columns:

- `DB_name`. The name of YOUR_PROJECT.
- `Sample`. Sample name of YOUR_METAGENOMITC_DATA.
- `CRISPR`. CRISPR id.
- `DRrepeat`.
- `Remark`. Upstream or downstrem of the CRISPR array.
- `Strand`.
- `Rank`. Distance bewteen the protein and the CRISPR array.
- `Prolen`. Protein size (aa).
- `Protein`. Protein id.
- `HMM`. HMMER annotation of this protein.
- `typetag`. Types of CRISPR-Cas systems identifed by cctyper (e.g. `I-*`, `II-*`,...). `Unkown` represents no CRISPR-Cas systems with certainty identified adjacent to the CRISPR array.
- `level`. Additional annotation of the protein:
  - `certain`: this protein was located within or very close to a known type of CRISPR-Cas system.
  - `putative`: this protein was located within or very close to a putative CRISPR-Cas system.
  - `large`: the size of this protein is over 2000aa.
  - `split`: this protein might be incomplete.
  - `other`: this protein was not included in the above mentioned scenarios.
  - `primary`: this protein falls under the 'other' category and has a size over 700 amino acids, indicating its potential to be a suspicious novel Cas protein.
- `pro_seq`. Amino acid sequence of the protein.
- `gene_seq`. CDS of the protein.
