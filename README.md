# CRISPRCasMiner

A CRISPR-Cas systems mining tools pipeline.

## Conda create environment
conda env create -f environment.yml
conda activate ccminer

## Run
chmod +x Path_to_CRISPRCasMiner/CRISPRCasMiner/ccminer/ccminer
Path_to_CRISPRCasMiner/CRISPRCasMiner/ccminer/ccminer -db_name Dataset_name --db DB_path --keep_tmp --cctyper_path path_of_cctyper_out Path_to_your_contigs output__path
 
