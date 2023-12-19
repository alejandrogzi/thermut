# thermut

This repository stores the code to reproduce the analysis of the project "The role of the thermodynamic pressure in Lensky's long-term evolution experiment"

To start this project clone this repo and do:

```bash
git clone https://github.com/alejandrogzi/thermut && cd thermut
./grouper.sh
```

to create the environment and build the population groups.

Your directory should now look like this:

```bash
...
├── groups
│   ├── m1.txt
│   ├── m2.txt
│   ├── m3.txt
│   ├── m4.txt
│   ├── m5.txt
│   ├── m6.txt
│   ├── m7.txt
│   ├── m8.txt
│   ├── p1.txt
│   ├── p2.txt
│   ├── p3.txt
│   ├── p4.txt
│   ├── p5.txt
│   ├── p6.txt
│   
...
```
Next, create the enviroment:

```bash
conda env create -f env.yml
conda activate thermut
```

Then, you are ready to run the pipeline:

```bash
nextflow run main.nf --dir /path/to/dir --meta /path/to/meta.csv --group /path/to/group.csv --gbk /path/to/gbk --fasta /path/to/fasta
```

Where:
- dir = directory containing fastq files
- meta = metadata file (./meta/metadata.txt)
- group = group file (a list of samples per population produced by ./grouper.sh)
- gbk = reference genbank file (./supp/REL606.6.gbk)
- fasta = reference fasta file (./supp/REL606.6.fasta)

