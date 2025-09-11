

# Overview

This repository contains the scripts used to analyze data and some data for the manuscript "Decoding Coral Resistance to Eutrophication: Selection of Hyper-Efficient Denitrifiers as Key Microbial Allies"

This repository appears to focus on **phylogenetic analysis, ancestral reconstruction, and gene arrangement studies** for microbial or genomic datasets.

## Main Scripts

- **ancestral_reconstruction.py**  

    This script performs large-scale **ancestral reconstruction** of orthologous genes across Ruegeria genomes, using both **AnGST** and **ecceTERA** reconciliation methods. The main workflow includes:

    1. **Phylogenetic tree preparation**  
    - Loads a species tree and prunes it to target genomes.
    - Selects representative genomes from each monophyletic cluster (based on nanopore sequencing, target list, or random selection).

    2. **Orthogroup inference**  
    - Runs **OrthoFinder** to define orthologous groups (OGs) from protein sequences.
    - Aligns OG sequences with **MAFFT** and builds gene trees with **IQ-TREE**.

    1. **Ancestral reconstruction**  
    - Runs **AnGST** and **ecceTERA** on each OG to reconstruct ancestral presence/absence and duplication/loss events.
    - Generates Slurm job arrays for HPC execution.

    1. **Post-processing and results**  
    - Compares OGs inferred by both methods, identifying shared and unique sets.
    - Annotates OGs with KEGG, InterPro, PANTHER, and CDD functional categories.
    - Outputs annotation tables (`ancestral_reconstruction/KEGG_OG_functions.tsv`).
    - Reconciliated Tree can be found in the folder `ancestral_reconstruction` named (`Angst_GeneTrees.tar.gz` and `ecceTERA_GeneTrees.tar.gz`)
    - used OG are here `ancestral_reconstruction/union_OG.faa`

    **Dependencies**: python packges: (`ete3`, `pandas`, `BioPython`, `tqdm`), OrthoFinder, MAFFT, IQ-TREE, AnGST, ecceTERA, InterProScan.


- **epddata.py**  
  This script processes **Environmental Protection Department (EPD)** water quality data, calculates mean nutrient concentrations for selected sampling stations, and converts units from mg/L to mM.

  - **Station mapping**: Maps original EPD station codes to project-specific codes.
  - **Data filtering**:  
    - Keeps only selected stations of interest.  
    - Filters by date range (≥ May 2010).  
    - Removes rows with missing values.
  - **Data selection**: Retains key water quality parameters (e.g., nitrogen species, phosphorus, turbidity, salinity, pH).
  - **Unit conversion**: Converts mg/L values to mM for each nutrient type.
  - **N:P ratio calculation**: Adds molar nitrogen-to-phosphorus ratio to the dataset.
  - **Output**: Writes results to an Excel file (`Table S1(mMversion).xlsx`).

- **genearrangment.py**  
  This script extracts and visualizes **denitrification gene clusters** from selected Ruegeria genomes, compares their gene arrangements, and highlights conserved genes across multiple genomes.

  - **Target region extraction**:  
  - Defines key denitrification genes (e.g., `nosZ`, `napA`, `norB`, `nirS`) from a reference genome.  
  - Extracts ~10 kb flanking regions from each target genome.
  - Saves extracted regions in FASTA (`.fna`), GenBank (`.gbk`), and protein FASTA (`.faa`) formats.
  - **Homology mapping**:  
  - Runs `blastp` between reference and target genome proteins to identify homologous genes.
  - Maps gene IDs to functional names based on the reference annotation.
  - **Synteny visualization**:  
  - Loads GenBank files into `pyGenomeViz` and draws gene arrows.
  - Highlights homologous denitrification genes in red.
  - Runs `blastn` between regions to detect conserved nucleotide segments and draws cross-links colored by identity.

- **genetrees.py**  
  This script automates **gene tree construction** for key denitrification genes and their auxiliary genes in Ruegeria genomes.  
  
  It extracts gene clusters from annotated genomes, filters for completeness, aligns sequences, builds phylogenetic trees, and generates iTOL annotation files.




---

## Folders

- **ancestral_reconstruction/**  
  Contains input data, results related to ancestral state reconstruction analyses.

- **figures/**  
  Stores generated figures (tree visualizations, heatmaps, synteny plots, etc.).

- **Trees_amplicon/**  
  Contains phylogenetic trees, built from amplifying region from all used Ruegeria genomes.

- **Trees_Dentrification/**  
  Contains phylogenetic trees using denitrification gene families from all used Ruegeria genomes.

- **phylosig_itol**  
  The results of phylogenetic signal (e.g., Blomberg's K, Pagel's λ) and prepare datasets for visualization in iTOL (Interactive Tree Of Life).
---

# Notes
Within the scripts, if you find any import like `from bin.format_newick import renamed_tree`. Please referred to the other repo `https://github.com/444thLiao/evol_tk`.


# Publication

Under review and submission.

# Data availability

The raw data for nirS amplicon and Ruegeria population resolving ATP5B, parC, and nirS amplicons are available under NCBI BioProject ID: PRJNA1310737, PRJNA1275610, PRJNA1275576. 

Raw reads and assembly of 419 Ruegeria genomes (i.e., missing raw reads for 10 genomes) were available under the NCBI BioProject ID [PRJNA1264799](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1264799?reviewer=mp7rulknapinfb95c9733nf72r). 

Raw reads and assembly of 34 Nanopore closed Ruegeria genomes were available under the NCBI BioProject ID [PRJNA1275854](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1275854?reviewer=g4mffnb8b965tq6b2662cpt6i8).


# Contact Us
If you have any questions or suggestions on these scripts, you are welcome to contact us via email: l0404th@gmail.com.
If you have any questions on the paper and experimental parts, you are welcome to contact us via email: hluo2006@gmail.com.