#!/bin/bash
######################################## Download DLPFC data ###################################################
## Download all preprocessed DLPFC slices
for fn in 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676; do
    wget -O ../data/DLPFC/${fn}_preprocessed.h5 https://zenodo.org/record/6334774/files/${fn}_preprocessed.h5?download=1
done

## Download DLPFC slice metadata
wget -O ../data/DLPFC/spatialLIBD_table_S5.csv https://zenodo.org/record/6334774/files/spatialLIBD_table_S5.csv?download=1

## Download raw sample III, slice 151674 for visualization-dlpfc.ipynb
download_slice_files() {
    wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/10X/$1/scalefactors_json.json -P sample-$2/$1/spatial
    wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/10X/$1/tissue_hires_image.png -P sample-$2/$1/spatial
    wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/10X/$1/tissue_lowres_image.png -P sample-$2/$1/spatial
    wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/10X/$1/tissue_positions_list.txt -P sample-$2/$1/spatial
    wget https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/$1_filtered_feature_bc_matrix.h5 -P sample-$2/$1
}
rm -rf ../data/DLPFC/sample-3
download_slice_files 151674 3
mv sample-3 ../data/DLPFC/
cd ../data/DLPFC/sample-3/151674/spatial
mv tissue_positions_list.txt tissue_positions_list.csv

## Download spot cell counts from original source
for fn in 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676; do
    wget --no-check-certificate --content-disposition -O ../data/DLPFC/${fn}_tissue_spot_counts.csv https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/Analysis/Histology/${fn}/tissue_spot_counts.csv
done

################################################################################################################

############################################# Download SCC data ###########################################
for n in 2 5 9 10; do
    wget -O ../data/SCC/scc_p${n}.zip https://zenodo.org/record/6334774/files/scc_p${n}.zip?download=1
    unzip -d ../data/SCC/ ../data/SCC/scc_p${n}.zip
done
rm -rf ../data/SCC/__MACOSX
rm ../data/SCC/*.zip

wget -O ../data/SCC/scc_visium.zip https://zenodo.org/record/6334774/files/scc_visium.zip?download=1
unzip -d ../data/SCC/ ../data/SCC/scc_visium.zip
rm -rf ../data/SCC/__MACOSX

###############################################################################################################

############################################# Download Stahl BC data ###########################################
wget -O ../data/Stahl-BC/stahl_bc_data.zip https://zenodo.org/record/6334774/files/stahl_bc_data.zip?download=1
## and unzip it
unzip -d ../data/Stahl-BC/ ../data/Stahl-BC/stahl_bc_data.zip
###############################################################################################################

############################################# Download Her2 data ###########################################
wget -O ../data/HER2/her2bc-ST-cnts.zip https://zenodo.org/record/6334774/files/her2bc-ST-cnts.zip?download=1
unzip ../data/HER2/her2bc-ST-cnts.zip -d ../data/HER2/
rm -rf ../data/HER2/__MACOSX
rm ../data/HER2/her2bc-ST-cnts.zip
mv ../data/HER2/ST-cnts/* ../data/HER2/
rm -rf ../data/HER2/ST-cnts
###############################################################################################################

############################################# Download Spinal-Cord data ###########################################
wget -O ../data/Spinal-Cord/spinal_cord_data.zip https://zenodo.org/record/6334774/files/spinal_cord_data.zip?download=1
unzip ../data/Spinal-Cord/spinal_cord_data.zip -d ../data/Spinal-Cord/
rm -rf ../data/Spinal-Cord/__MACOSX
rm ../data/Spinal-Cord/spinal_cord_data.zip
###############################################################################################################


############################# Download and Process Misc. like caches ###################################
## Unzip included stutility alignment cache
unzip -d ../data/DLPFC/saved_results/ ../data/DLPFC/saved_results/stutil-alginments.zip

wget -O ../data/DLPFC/saved_results/tangram_pis.pickle https://zenodo.org/record/6395124/files/tangram_pis.pickle?download=1

wget -O ../data/DLPFC/saved_results/center2_a0.1_KL_seed2.h5ad https://zenodo.org/record/6395124/files/center2_a0.1_KL_seed2.h5ad?download=1

wget -O ../data/DLPFC/saved_results/sample3-X-all-seurat-integrated.h5ad https://zenodo.org/record/6395124/files/sample3-X-all-seurat-integrated.h5ad?download=1
