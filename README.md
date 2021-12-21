# siRNA
# This code is a used for data processing and analysis of the siRNA imaging data produced by the microscope.
# Data generated in primary screen with Whole genome Ambion library plates read by ScanR microscope was combined with the annotations of library plates. As a part of quality control check, distribution and performance of the controls was analysed per plate with the help of the well signals from wells marked with positive, negative controls and internal controls. For each of the plate, ratio of the baku_green_objects to Dapi_Objects was calculated. To take care of plate wise differences, Percent to Negative Control normalization was used (Birmingham et al,) in which the bac_nuc_ratio values were normalized to the median value of wells marked with the scrambled siRNA controls (negative control). The resulting ratios were used to filter and select interesting hits for the secondary screen.
# In the secondary screen along with the negative controls, the bac_nuc_ratio data was normalized to the median value of the wells which had been added Antibodies (Antibody block).
# The data was processed with the help of R programming language and the plots were generated using GGPlot2 R package. 

# Reference : Birmingham, Amanda et al. “Statistical methods for analysis of high-throughput RNA interference screens.” Nature methods vol. 6,8 (2009): 569-75. doi:10.1038/nmeth.1351
