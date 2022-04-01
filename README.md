# Spreading_Lensch_2022
Code and data associated with Dynamic spreading of chromatin-mediated gene silencing and reactivation between neighboring genes in single cells: 
https://www.biorxiv.org/content/10.1101/2021.11.04.467237v1.full

## CUT&RUN Analysis
Raw sequencing data from CUT&RUN experiments can be found on GEO server: GSE189540. 

Accessory scripts "edit_sam_*.py" are edited to remove ambiguous reads that fall entirely within duplicated regions 

Fastq files are processed first with "cutrun_pipeline_v9_2021-10-16.py" to align, filter and normalize reads, followed by "cutrun_pipeline_part2_20211019.py" to generate bedgraph files. 

Shell scripts "find_5kbreporter_lines.sh" or "find_chr19_lines.sh" for CHO or K562, respectively, are used to generate smaller bedgraphs for generating quantification plots. 

Python scripts "bedgraphAUC_K562_5kb.py", "bedgraphAUC_K562_5kb_signal_diff_norm.py", "catplot_K562_5kb_signal_diff_norm" and "barplot_K562_5kb_signal_diff_norm.py" are used to generate various plots quanitifying signal across elements of interest.  

## Time-lapse microscopy analysis
Single cells in acquired TIFF images are segmented and tracked MACKtrack (https://github.com/brookstaylorjr/MACKtrack). 

Matlab script "cell_analysis_macktrack_v3.0.m" is then used to filter pre-silenced cells or poorly tracked cell traces. This script is dependent on "copychildren.m" in the same directory. 

Matlab script "trace_stitchin_script_v3.0.m" is then used to call cell division event and stitch traces together in order to call silencing events.

## Model for Time Evolution of Fluorescence Distributions
Raw data exists as .csv files named for CR+Spacer/Geometry+Rep/Clone (K562 reps are distinguished by an "R" before the rep number).
In each .csv file, there are 3 columns for the following fields: [log10(Cit_Fluorescence), log10(mCh_Fluorescence), Day].
Data was gated for live cells and IFP pos using EasyflowCyto.m before writing to .csv format. No dox days are negative integers.

Fits are conducted in "Probabilistic Model Fitting..." files. All fits include mRNA and spike delay, but HDAC4 fits have an unfixed
gamma while KRAB fits use gamma values from HDAC4 fits, or estimated gamma values from raw HDAC data. The code for each core file
is identical, but the fits are distributed among the core files to allow parallel processing.

The best fit parameters are saved as .pkl files, which can be opened in Python.
The returned parameters are listed in the following order:
[alpha, beta, gamma, sigma, k, mu, P_b, F_b, sd]

Violin plots and parameter plots are generated in "Half Violin Plots + Model including mRNAsd-Final" and "Prob Model Param Plots with 90% CI Final" Jupyter notebooks, respectively.

