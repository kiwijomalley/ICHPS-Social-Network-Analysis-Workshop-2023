# ICHPS-Social-Network-Analysis-Workshop-2023

Workshop to be taught at the 2023 International Conference on Health Policy Statistics in Scottsdale Arizona, USA

The following files are supplemental to the workshop. Participants are welcome to review these and experiment with the code and data prior to the workshop. 

ICHPS_ShortCourse_2023: Powerpoint slides containing expanded descriptions and extracts of code and output for the examples discussed in the workshop.
WorkshopICHPS2023: R Code for performing descriptive analyses presented in workshop for the physician clinic network example
WorkshopICHPS2023models: R Code for estimating statistical models of the network relational data and of peer-associations (social influence or contagion) presented in workshop for the physician clinic network example
WorkshopICHPS2023grandpa: R Code for performing descriptive analyses and for estimating statistical models of relational data as well as peer associations for the physician network embedded in the grandpa package
ICHPS_ClinEdgeList.csv: Edgelist for the physician clinic network
ICHPS_PhysPract.txt: Adjacency matrix for the physician clinic network
ICHPS_nodecov.dat: Attribute data set for the physicians in the physician clinic network

The physician clinic network was described in Keating et al (2007) and re-analyzed in O'Malley and Marsden (2008). It is a small network containing only 33 actors. The network in the grandpa package is based on a real bipartite patient-physician encounter data set. However, for confidentiality, the actual data is replaced with a randomly drawn replica based on the multiple attribute augmented model described in Bobak et al (2022).

The R scripts each contain multiple analyses, allowing users to extract and adapt segments of code for their own purposes. They are heavily commented and should be self explanatory to follow, especially with the aid of the slides from the short course.

References:
Keating, N. L., Ayanian, J. Z., Cleary, P. D., and Marsden, P. V. (2007), Factors affecting influential discussions among physicians: a social network analysis of a primary care practice," Journal of general internal medicine, 22, 794-798
O'Malley, A. J. and Marsden, P. V. (2008), The Analysis of Social Networks, Health Services & Outcomes Research Methodology, 8, 222-269
Bobak CA, Zhao Y, Levy JJ, and O’Malley AJ. (2022). GRANDPA: GeneRAtive Network sampling using Degree and Property Augmentation applied to the analysis of partially confidential healthcare networks. arXiv.2211.15000
