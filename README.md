# Visualize Nigon units in Nematode assemblies

You need to run BUSCO on your genome assembly using the nematoda_odb10 dataset.
This script takes as input the resulting `full_table.tsv`.


## Installation
The script relies on several R libraries. These can be installed with conda. I recommend creating a new conda environment to avoid conflicts with other packages.

```
# Define name for new conda environment
ENV_NAME=vis_alg
# Create environment and install dependencies
conda create -n $ENV_NAME -c r -c conda-forge icu=58 r-dplyr=0.8 r-readr=1.3.1 r-scales=1.1.1 r-gtools=3.8.2  r-optparse=1.6.6 r-ggtext=0.1.0 xorg-libxrender=0.9.10 -y
# Activate conda environment
conda activate $ENV_NAME
```



## Usage
There are two examples which you can use to test.

```
git clone https://github.com/pgonzale60/vis_ALG.git
cd vis_ALG/

Rscript bin/vis_ALGs.R -b full_table.tsv -n data/gene2Nigon_busco20200927.tsv.gz -s Genus_species -o output.png
```

## Examples

*Caenorhabditis elegans* and *Oscheius tipulae* full_table.tsv resulting from BUSCO 4 using nematoda_odb10 are included as examples.

To generate a PNG image of Nigon units in *O. tipulae* chromosomes you can execute
```
Rscript bin/vis_ALGs.R -b examples/oscheius_tipulae.local.v3_1rx.full_table.tsv -n data/gene2Nigon_busco20200927.tsv.gz -s Oscheius_tipulae -o otipu.png
```

To generate a JPEG image of Nigon units in *C. elegans* chromosomes you can execute
```
Rscript bin/vis_ALGs.R -b examples/oscheius_tipulae.local.v3_1rx.full_table.tsv -n data/gene2Nigon_busco20200927.tsv.gz -s Caenorhabditis_elegans -o cele.jpeg
```
