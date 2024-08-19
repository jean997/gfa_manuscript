#!/bin/bash

mkdir -p gwas_data

wget https://zenodo.org/records/13250706/files/ld_reference_data.tar.gz 
tar -xzf ld_reference_data.tar.gz
rm ld_reference_data.tar.gz

wget https://zenodo.org/records/13250706/files/autoimmune.tar.gz 
wget https://zenodo.org/records/13250706/files/bloodcell.tar.gz 
wget https://zenodo.org/records/13250706/files/t2d_cad.tar.gz 


tar -xzf autoimmune.tar.gz
rm autoimmune.tar.gz
mv autoimmune gwas_data/

tar -xzf bloodcell.tar.gz
rm bloodcell.tar.gz
mv bloodcell gwas_data/

tar -xzf t2d_cad.tar.gz
rm t2d_cad.tar.gz
mv t2d_cad gwas_data/
