# ICHPS-Social-Network-Analysis-Workshop-2023

Workshop on Social Network Analysis for 2023 International Conference on Health Policy Statistics in Scottsdale Arizona, USA. Presented January 9, 2023 by James O'Malley (James.OMalley@Dartmouth.edu)

The following files are supplemental files for the workshop. Participants are welcome to review these and experiment with the code and data prior to the workshop. 

ICHPS_ShortCourse_2023.pdf: Slides containing expanded descriptions including extracts of code and output for the examples discussed in the workshop

WorkshopICHPS2023.R: R Code for performing descriptive analyses presented in workshop for the physician clinic network example

WorkshopICHPS2023models.R: R Code for estimating statistical models of the network itself and of models of peer-associations (social influence or contagion) presented in workshop for the physician clinic network example

WorkshopICHPS2023grandpa.R: R Code for performing descriptive analyses, for estimating statistical models of the network itself, and for modeling peer associations for the physician network embedded in the grandpa package of Bobak et al (2022)

ICHPS_ClinEdgeList.csv: Edgelist for the binary-valued directed physician clinic network

ICHPS_PhysPractBin.txt: Adjacency matrix for the physician clinic network

ICHPS_nodecov.dat: Physician attribute data set for the physician clinic network

The physician clinic network is described in Keating et al (2007) and re-analyzed in O'Malley and Marsden (2008). It is a small network containing 33 physicians. The network in the grandpa package is based on a bipartite patient-physician encounter network data set. However, for confidentiality, the actual network (and attributes) are replaced with an approximate randomly drawn network of the same size generated using the augmented multiple attribute network model described in Bobak et al (2022) and presented in a poster at the conference.

The R scripts contain multiple analyses. Users may extract and adapt segments of code for their own purposes. They are heavily commented and should be self explanatory to follow, especially with the aid of the slides from the workshop.

References:

Keating, N. L., Ayanian, J. Z., Cleary, P. D., and Marsden, P. V. (2007), Factors affecting influential discussions among physicians: a social network analysis of a primary care practice, Journal of general internal medicine, 22, 794-798

O'Malley, A. J. and Marsden, P. V. (2008), The Analysis of Social Networks, Health Services & Outcomes Research Methodology, 8, 222-269

Bobak, C. A., Zhao, Y., Levy, J. J., and O’Malley, A. J. (2022). GRANDPA: GeneRAtive Network sampling using Degree and Property Augmentation applied to the analysis of partially confidential healthcare networks. arXiv.2211.15000
