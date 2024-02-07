#!/usr/bin/env python
import pandas as pd
import scanpy as sc
import numpy as np
import os
import pegasusio as io
import argparse

parser = argparse.ArgumentParser(description='Parser for Downsampling - Demultiplexing')
parser.add_argument('--hto_data', help='path to HTO data')
parser.add_argument('--downsampling_seed', help='Number used as seed for the process',type=int, default=32)
parser.add_argument('--downsampling_percentage', help='Percentage used for downsampling as int number. Ex:30',type=int, default=30)
parser.add_argument('--hto_data_name', help='name to save downsampled HTO data')
parser.add_argument('--path_save_hto', help='path to save HTO  dowsampled data')
parser.add_argument('--path_hto_h5', help='path to save HTO  dowsampled data H5')

args = parser.parse_args()

def downsample_input(hto_path,fix_seed, percentage):

    print("seed is")
    print(fix_seed)
    print("percentage is")
    print(percentage)
    if os.path.isdir(hto_path):
		print(f"{hto_path} mtx HTO matrix provided")
		hto_data = sc.read_10x_mtx(hto_path,gex_only=False)
    elif os.path.isfile(path):
		print(f"{hto_path} HTO H5 object")
		hto_data = sc.read_10x_h5(hto_path,gex_only=False)
    sampling_percentage =  percentage
    seed = fix_seed
    np.random.seed(seed)
    num_rows_to_sample = int(hto_data.shape[0] * (sampling_percentage / 100))
    print(f"rows sampled: {num_rows_to_sample}")
   # Sample rows
    sampled_indices = np.random.choice(hto_data.shape[0], num_rows_to_sample, replace=False)
	sampled_hto = hto_data[sampled_indices, :]
    return sampled_hto



percentage = args.downsampling_percentage
set_seed = args.downsampling_seed

##Downsampling
sampled_hto = downsample_input(args.hto_data,set_seed,percentage)

##The creation for intermediate files is necessary to create compatible files for Hadge and its inpput tools
##Intermediate files with downsampled data creation

path_save_hto =  args.path_save_hto

#write  H5ad Anndata object
sampled_hto.write_h5ad(path_save_hto)

#path to save H5
path_hto_h5  = args.path_hto_h5

#read anndata sampled for h5
hto_h5 = io.read_input(path_save_hto)

#write h5 
io.write_output(hto_h5,output_file = path_hto_h5, file_type="h5")



