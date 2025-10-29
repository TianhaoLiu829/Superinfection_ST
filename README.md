# Influenza, Staphylococcus aureus Super-infection Disrupts Spatially Coordinated Cellular Immunity in the Lung
# abstract
Influenza-associated bacterial super-infections in the lung lead to increased morbidity and mortality. Previous studies have investigated how preceding viral infection causes dysregulation of the innate and adaptive immune systems, leading to increased susceptibility of developing secondary bacterial pneumonia. However, these previous studies cannot account for the spatial context of immune cells in lung. In our study, we employed a spatial transcriptomics platform (10X Genomics Visium) to systematically characterize the coordination of immune cells during super-infection. We compared deconvoluted spatial transcriptomics data between super-infection and single influenza and methicillin-resistant Staphylococcus aureus (S. aureus) infections. Consequently, we found that the recruitment of neutrophils and interstitial macrophages from lung parenchyma to the airways was inhibited in super-infection, likely impairing pathogen clearance. Additionally, by analyzing cell colocalization and signaling, we found that the interaction between CD4+ T cells, B cells, and dendritic cells was disrupted by secondary bacterial super-infection. Our study constructs a spatial sequencing atlas of lung super-infection, highlighting how secondary bacterial challenge significantly impacts recruitment and signaling of immune cells. These findings provide insight into the aberrant inflammation in the super-infected lung and may aid development of therapeutics that target key immune cell recruitment pathways to improve host response to super-infection.
# Code for analysis
The code was splited into three sessions (original_analysis, figure and object).

The folder of "original_analysis" contains code for (1) RCTD cell type deconvolution, (2) sptial communication, (3) Unsupervised clustering and (4) annotation of airway and inflammatory regions. The code in the folder of original_analysis should be ran in the order above. Notably, some spots fail in RCTD or commot will be first excluded from the study.

The folder of "object" show be dowloaded from zenodo. The files contain the result from original analysis. The objects can be directly used to create figures with the code in folder of "figure".

The folder of "figure" contains code for making main figures and supplementory figures. 
<img width="3994" height="2181" alt="1 - Figure1" src="https://github.com/user-attachments/assets/980e0f44-2b6b-444b-b7c8-c2674ca58ec7" />

