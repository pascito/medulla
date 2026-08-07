import os
import sys

# Setup Server Paths
spineplot_root = '/exp/icarus/app/users/aanilkum/medulla/spineplot'
os.chdir(spineplot_root)
sys.path.insert(0, spineplot_root) 
from analysis import Analysis

# Point to your single, updated TOML and ROOT file
data_file = '/exp/icarus/data/users/aanilkum/medulla/output_root/1muNp0pi_Nge1_uncontained_batch_dpT_sys.root'
main_toml = '/exp/icarus/app/users/aanilkum/medulla/spineplot/My_analysis/dPT/dpT.toml'

print("--- Generating ALL masked plots with full systematics ---")
my_analysis = Analysis(main_toml, data_file)
my_analysis.run(close_figs=False)
print("--- Finished! ---")