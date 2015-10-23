

########################################################################################################################################
################################################## How to run GOAna ####################################################################
########################################################################################################################################

#+++++++++++++++++++++++++++++++++++++++++++++++++++ Preparation ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

Packages which need to be installed:
- affy
- annotate
- GO.db
- annotation library for the Array Chip in use (e.g. hgu133plus2.db or hgu133a.db and so on)

#+++++++++++++++++++++++++++++++++++++++++++++++++++ Description +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

Provided scripts:
- defineGOIDs.R
- getdistances_calls.R
- CreateCytoscape.R

defineGOIDs:
- defines the list of GOIDs the analysis will be run on
- Input: - chip:               Array type (without ".db" at the end)		default: hgu133a
         - numprobesets_low:   the minimal number of probesets per GOID		default: 5
         - numprobesets_high:  the maximum number of probesets per GOID		default: 10000
         - category:           GO category of IDs							default: "BP" (biological processes)
         - GO.node:            point in the GO hierarchy to start			default: "GO:0008150" (biological_process)
 - Output: vector of GO IDs together with corresponding probe IDs

getdistances_calls:
- calculates distances 
- Input:  - inputdata: 	expression matrix with probesets as rownames
	      - calls: 		T or F (TRUE or FALSE), whether there is a call (A,P,M) for each ProbeID on the chip
		  - call.file: 	if calls == T, then you have to give the file which includes the calls						default: TRUE
		  - numpresent: minimum number of samples in which a probe should be present								default: 5
          - GOIDs: 		GOIDs to use as returned by defineGOIDs
          - classes: 		sample affiliations
          - outputfile: 	name of output file
          - chip:			array type to use, e.g. "hgu133a", "illuminaHumanv3"									default: "hgu133a"
          - permutations: number of permutations to calculate significance											default: 1000
- Output: .txt file containing the GOIDs together with p-values; in the case of more than two groups,
           the file will contain a p-value for each pairwise comparison
 
 
CreateCytoscape
- extends the significant GOIDs and creates a contrib-network 
- Input: - GOIDs: a vector of GOIDs to analyze
 	     - calls: T or F (TRUE or FALSE), whether there is a call (A,P,M) for each ProbeID on the chip		default: TRUE
         - datacalls: if calls == T, then you have to give the file which includes the calls
         - intersection: determines the number of shared GOIDs between two genes							default: 2
         - chip: Array type to use, e.g. "hgu133a", "illuminaHumanv3"										default: "hgu133a"
		 - comparison: in the case of more than two groups, specify the pairwise comparison of interest 	
						as character, it will be used for the output file names	(e.g. "A_vs_B")				default: "1"
 
 
 function(GOIDs,calls=T,datacalls,intersection=2,chip=, )
 

#++++++++++++++++++++++++++++++++++++++++++++++++++++ Example ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# set working directory
setwd("X/Y/Z")

# load script source files
source("defineGOIDs.R")
source("getdistances_calls.R")
source("CreateCytoscape.R")

# load example data
# reference: Mootha et al. 2003, Nature Genetics)
load("diabetes_data.RData")
load("classes_diabetes_data.RData")

# calculate GOIDs
GOIDs_BP <- defineGOIDs(chip = "hgu133a")

# calculate distances, this example runs about (start: 11.28)
distances <- getdistances_calls(diabetes_data,calls=F,GOIDs=GOIDs_BP,classes=classes_diabetes_data,outputfile="ResultGOAna_DiabetesData_BP.txt", chip = "hgu133a")

# read in results text file
res_GOAna <- read.delim("ResultGOAna_DiabetesData_BP.txt", header = TRUE, stringsAsFactor = FALSE, check.names = FALSE)

# Define the column name for the p-values of the comparison you are interested in
group <- "p-value 1 vs 2"

# define significant GO IDs with a desired p-value cutoff, here 0.01
sig_GOIDs <- res_GOAna[res_GOAna[,group] < 0.01,1]

# create a Cytoscape network 
createCytoscape(sig_GOIDs,calls=F, chip = "hgu133a", comparison = group)
# depending on the data, the network will be quite large; increase the interaction parameter to reduce the network

# Finally, import the generated .sif file as network into Cytoscape



