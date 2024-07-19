Simple introduction to files and folders.

Packages "ggplot2", "gridExtra", "glmnet" are used.
Sources install_packages.R to check whether they are installed. If not, the program will automatically install the packages.

algorithm_on_maximin.R, COA_construct.R, functions.R are three core files.
-------------------------------------------------------------------

algorithm_on_maximin.R: provide algorithm to search for MCOA.

COA_construct.R: construct initial COA.

functions.R: transform design to location matrix and model matrix.

-------------------------------------------------------------------



For convenience, we have copied the three important files into each individual leave folder.
----------------------------------------------------------------------------------------

maximin design: generate maximin distance design in our paper. See README in the folder.

MCOA: all MCOAs in our paper.

MCOAl2: MCOAs obtained under the L2-distance.

COAD: COA with high PWO D-efficiency.

workcost: simulation of Case study. See README in the folder.

shapley: simulation of Application to Shapley value. See README in the folder.

----------------------------------------------------------------------------------------