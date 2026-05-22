
import scanpy as sc
import numpy as np
import loompy as lp
import sys
import os

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <input.csv> <output.loom>")
    sys.exit(1)

input_csv = sys.argv[1]
output_loom = sys.argv[2]

output_dir = os.path.dirname(output_loom)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir)

# read csv file
adata = sc.read_csv(input_csv)
print("Loaded data:")
print(adata)

# create loom file
row_attrs = {"Gene": np.array(adata.var_names)}
print(row_attrs)
col_attrs = {"CellID": np.array(adata.obs_names)}
print(col_attrs)
lp.create(output_loom, adata.X.transpose(), row_attrs, col_attrs)
print(f"Successfully saved loom file to: {output_loom}")

## Example usage:
# python csv2loom.py input.csv output.loom
