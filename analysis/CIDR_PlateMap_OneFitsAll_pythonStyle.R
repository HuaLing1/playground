#!/mnt/linuxtools/R/3.5.1/bin/Rscript

library(optparse)
library(dplyr)
library(xtable)
library(PlateMap)
options(stringsAsFactors = FALSE)
rm(list=objects())


##########
# 0. Parsing arguments 
##########
option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder name"),
  make_option(c("-i", "--inputpath"), type="character", default=NULL, 
              help="path to input files"),
  make_option(c("-a", "--array"), type="character", default=NULL, 
              help="input file name for rearray file"),
  make_option(c("-d", "--dup"), type="character", default=NULL, 
              help="input file name for duplicate pair file"),
  make_option(c("-f", "--fenotypeOrSIF"), type="character", default=NULL, 
              help="input file name for SIF file"),
  make_option(c("-n", "--numSamples"), type="integer", default=NULL, 
              help="number of samples"),
  make_option(c("-s", "--stratified"), type="character", default=NULL, 
              help="variable list for stratified varaibles"),
  make_option(c("-b", "--balanced"), type="character", default=NULL, 
              help="variable list for checking balance after randomization"),
  make_option(c("-r", "--random"), type="integer", default=NULL, 
              help="random seeds")
); 
length(option_list)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


cat("proj:\n")
cat(opt$p)
cat("\n\n")
cat("inputPath:\n")
cat(opt$i)
cat("\n\n")
cat("arrayFile:\n")
cat(opt$a)
cat("\n\n")
cat("dup:\n")
cat(opt$d)
cat("\n\n")
cat("fenotypeOrSIF:\n")
cat(opt$f)
cat("\n\n")
cat("numSamples:\n")
cat(opt$n)
cat("\n\n")
cat("stratified:\n")
cat(opt$s)
cat("\n\n")
cat("balanced:\n")
cat(opt$b)
cat("\n\n")
cat("random:\n")
cat(opt$r)
cat("\n\n")


cat("\n\n
###################################################################
#    1. Set up working directory and project specific variables   #
###################################################################
\n\n")
proj <- opt$p
cmd <- paste0("mkdir /mnt/research/statgen/UW_PlateMap/", proj); system(cmd)
cmd <- paste0("mkdir /mnt/research/statgen/UW_PlateMap/", proj, "/data"); system(cmd)
cmd <- paste0("mkdir /mnt/research/statgen/UW_PlateMap/", proj, "/output"); system(cmd)
dir.work <- file.path("/mnt/research/statgen/UW_PlateMap/", proj)

dir.input <- opt$i
file.rearray <- opt$a
file.dup <- opt$d
file.SIF <- opt$f

# type in number of samples
nsamp <- opt$n


# randomized variable list for check
# vars <- c("Sex", "Population", "DNA.Source", "DNA.Extraction.Method", "MTA", "IRB")
if(unlist(opt$s)=="NA"){
  vars <- unlist(strsplit(unlist(opt$b), " "))
} else if (unlist(opt$s)!="NA"){
  vars <- c(unlist(strsplit(unlist(opt$s), " ")), unlist(strsplit(unlist(opt$b), " ")))
}

randomNum <- opt$r



cat("\n\n
#############################################################################
#   2. Read in 3 input files: 1) Duplicate pair 2) SIF_Aliquot 3) Rearray   #
#############################################################################
\n\n")
input.dup <- read.csv(file.path(dir.input, file.dup), as.is=T, header=T); dim(input.dup)
# [1] 69   12
table(input.dup$Expected_Aliquot_1.Status, input.dup$Expected_Aliquot_2.Status)
#           Canceled Created
# Canceled        1       0
# Created         4      64
col.subjectID <- grep("SubjectID", names(input.dup))
head(input.dup)[col.subjectID]
# SubjectID_1 SubjectID_2
# 1       C2_49     C2_9902
# 2     C1_2947     C1_9970
# 3       D_302      D_9901
# 4       C2_19     C2_9901
# 5     C2_1848     C2_9941
# 6      C2_831     C2_9918
##### Comments: 
# 1) 69 Dup pairs, AND 64 were Active. 
# 2) Dup pairs have different subjectID for this project. 

temp1 <- read.csv(file.path(dir.input, file.SIF), as.is=T, header=T); 
input.SIF <- temp1[!temp1$DNA_Class %in% c("0", "null"),]
dim(input.SIF)
# [1] 3233   28
names(input.SIF)
# [1] "Family"                       "Individual"                   "Father"                      
# [4] "Mother"                       "Sex"                          "DNA_Class"                   
# [7] "IRB"                          "MTA"                          "Population"                  
# [10] "Organism"                     "DNA.Source"                   "DNA.Extraction.Method"       
# [13] "Subject_ID"                   "Expected.Aliquot"             "Expected.Aliquot.Status"     
# [16] "Aliquot.ID"                   "Aliquot_Status"               "SIF.Investigator.Column.1"   
# [19] "SIF.Investigator.Column.2"    "SIF.Investigator.Column.3"    "SIF.Investigator.Column.4"   
# [22] "SIF.Investigator.Column.5"    "MANREC.Investigator.Column.1" "MANREC.Investigator.Column.2"
# [25] "MANREC.Investigator.Column.3" "MANREC.Investigator.Column.4" "MANREC.Investigator.Column.5"
# [28] "X"

# constrcut Duplicate list. This depends on whether the duplicate pairs have same or different subject ID
Subject.ID.1 <- input.SIF$Subject_ID[duplicated(input.SIF$Subject_ID)]
(ct.dup <- length(Subject.ID.1)); (ct.dup.uniq <- length(unique(Subject.ID.1)))
# [1] 0 # this is to check whether all duplicates are doubletons. if ct.dup==ct.dup.uniq!=0, this means all dup are doubletons
# if ct.dup==0, then dup has diff subjectID. 
if (ct.dup==0){
  dups <- as.data.frame(cbind(input.dup$SubjectID_1[input.dup$Pair.Status=="Active"],
                              input.dup$SubjectID_2[input.dup$Pair.Status=="Active"]))
  input.SIF$SampleID <- input.SIF$Subject_ID
} else if (ct.dup==ct.dup.uniq) {
  Subject.ID.2 <- paste(Subject.ID.1, "_DUP", sep="")
  dups <- as.data.frame(cbind(Subject.ID.1, Subject.ID.2))
  input.SIF$SampleID <- ifelse(duplicated(input.SIF$Subject_ID), 
                               paste(input.SIF$Subject_ID, "_DUP", sep=""), 
                               input.SIF$Subject_ID)
} else {
  print("non-doubleton duplicates")
}
dim(dups) #[1] 64  2 ##### comments: This is same number as the active pair from dup file. 
write.csv(dups, file = file.path(dir.work, "data/dups_reformated.csv"), row.names = FALSE, quote = FALSE)

# create new variable called SampleID so that duplicates have different subject ID. 
table(duplicated(input.SIF$SampleID)) # Expect to see all False
# FALSE 
# 3233 



cat("\n\n
#########################
#   3. MAKE PLATEMAP    #
#########################
\n\n")
# define the strata
# stratified variable list
if(unlist(opt$s)=="NA" & !all(input.SIF$Family==0) ){
  var.strata <- "Family"
} else if (unlist(opt$s)=="NA" & all(input.SIF$Family==0)){
  var.strata <- NULL
} else if (unlist(opt$s)!="NA" & !all(input.SIF$Family==0)){
  var.strata <-c(unlist(strsplit(unlist(opt$s), " ")), "Family")
} else if (unlist(opt$s)!="NA" & all(input.SIF$Family==0)){
  var.strata <- unlist(strsplit(unlist(opt$s), " "))
}


# check if there is any cross-study dup. If yes, create Reserve for DNA_Class=3
if(any(input.SIF$DNA_Class==3)){
  input.SIF$Reserve <- input.SIF$DNA_Class==3
  sample.strata <- input.SIF[ , c("SampleID", var.strata, "Reserve"), drop=FALSE]; dim(sample.strata)
} else {
  sample.strata <- input.SIF[ , c("SampleID", var.strata), drop=FALSE]; dim(sample.strata)
}
input.rearray <- read.csv(file.path(dir.input, file.rearray), as.is=TRUE, header=TRUE, skip=5)
dim(input.rearray) 
# 3648   6
names(input.rearray)
# [1] "Rearray_Container_ID"                                      
# [2] "Well"                                                      
# [3] "Subject_ID..must.match.the.Subject_ID.provided.in.the.SIF."
# [4] "CIDR.Comment"                                              
# [5] "Investigator_Column_1"                                     
# [6] "Investigator_Column_2"   
table(input.rearray[,5])
table(input.rearray[,6])
# all NA

# change names to be the ones expected by the PlateMap function
names(input.rearray)[1:3] <- c("Plate", "Well", "SampleID")
# check the number of CIDR controls per plate
ct.control.plate <- as.data.frame(table(input.rearray$Plate, input.rearray$SampleID!="")) %>%
  filter(Var2 == T) 
lower.ct <- range(ct.control.plate$Freq)[1]; higher.ct <- range(ct.control.plate$Freq)[2]
stopifnot(lower.ct == higher.ct)
if( lower.ct == higher.ct ){
  num.control.plate <- lower.ct
  print(paste0("number of CIDR controls on each plate is: ", lower.ct))
} else {
  num.control.plate <- round(sum(ct.control.plate$Freq)/dim(ct.control.plate)[1], 2)
  cat(paste0("number of CIDR controls by plate ranges from ", lower.ct, " to ", higher.ct,
             "\nAvg control ct per plate is: ", num.control.plate))
}

# if there are too many samples, we sometimes have to remove the last plate
(n_plates_to_keep = ceiling(nrow(sample.strata) /(96-num.control.plate))) # 35
# figure out how many plates we should have
(n_plates_rearray <- length(unique(input.rearray$Plate))) # 38
# 2 CIDR controls on each plate
# remove the last several plates
n_plates_remove <- n_plates_rearray - n_plates_to_keep
plates <- unique(input.rearray$Plate)
lastplates <- rev(plates)[0:n_plates_remove]
input.rearray <- input.rearray[!(input.rearray$Plate %in% lastplates), ]
dim(input.rearray) # [1] 3360    6

# save header for later
con <- file(file.path(dir.input,file.rearray), "r")
cidrHeader <- readLines(con, n=6)
close(con)

# Run plate map function
set.seed(randomNum) # run it everytime you call plateMap() function if you want to get the same map file
pm <- plateMap(sample.strata, input.rearray, duplicates=dups)
dim(pm)
# [1] 3360    6

# write out new plate map
pm[is.na(pm)] <- ""
# wirite one for checking purpose
pm.filename <- paste0("plateMap_for_checking_", Sys.Date(), ".csv")
pm.file <- file.path(dir.work, "data", pm.filename)
con <- file(pm.file, "w")
writeLines(cidrHeader, con)
write.table(pm, con, sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
close(con)

# remove "_DUP" at the end, and used for produciton
SampleID <- gsub("_DUP", "", pm$SampleID)
pm$SampleID <- SampleID
pm.filename2 <- paste0(proj, "_RearrayManifest_", Sys.Date(), ".csv")
pm.file2 <- file.path(dir.work, "output", pm.filename2)
con <- file(pm.file2, "w")
writeLines(cidrHeader, con)
write.table(pm, con, sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
close(con)



cat("\n\n
#########################
#   4. CHECK PLATEMAP   #
#########################
\n\n")

##### (1) Original manifest check
pm.file <- file.path(dir.work, "data", pm.filename)
pm <- read.csv(pm.file, as.is=TRUE, skip=5)
dim(pm); names(pm); head(pm) # [1] 11904     6
names(pm)[3] <- "Subject_ID"
plate.manifest <- read.csv(file.path(dir.input, file.rearray), as.is=TRUE, header=TRUE, skip=5)
# add some code here to handle check of additional plates
plate.manifest.removed <- plate.manifest[! plate.manifest$Rearray_Container_ID %in% lastplates,]
stopifnot(all.equal(pm[,c(1,2,4)], plate.manifest.removed[,c(1,2,4)])) # T
ctrl <- grep("CONTROL", pm$Subject_ID)
stopifnot(all.equal(pm$Subject_ID[ctrl], plate.manifest$Subject_ID[ctrl])) # T
# new plate map matches original CIDR manifest

# check plates and wells
names(pm)[1] <- "Plate"
table(pm$Plate)
table(table(pm$Plate)) # 35 batches of 96 wells
table(pm$Well)
table(table(pm$Well)) # 96 wells on 35 plates

# check for NA or empty cells
tmp <- unlist(lapply(pm, function(x) sum(is.na(x)))); tmp[tmp!=0]
# Investigator columns
tmp <- unlist(lapply(pm, function(x) sum(x %in% ""))); tmp[tmp!=0]
# Subject_ID CIDR.Comment 
# 57        3290 

# check that numbers of wells add up
(nwells <- nrow(pm)) # 6816
(nctrl <- length(ctrl)) # 142

# empty wells
idxEmpty <- which(pm$Subject_ID == "")
(nempty <- length(idxEmpty)) # 61
stopifnot(nsamp + nctrl + nempty == nwells)
# are empty wells all at the end?
range(idxEmpty) # [1] 6755 6816
cat("##### Manifest check is done \n\n")


##### (2) phenotype check 
# merge the plate map with the phenotype data
names(input.SIF)[13] <- "Subject_ID"
pm.phen <- pm %>% 
  inner_join(input.SIF, by = "Subject_ID")

# check chi-sq values here
# but print the tables in the Rmd
# vars <- c("Population", "DNA.Source", "DNA.Extraction.Method", "MTA", "IRB") # this was moved to the beginning of the script
for (var in vars) {
  x <- table(pm.phen$Plate, pm.phen[[var]])
  pval_all <- suppressWarnings(chisq.test(x)$p.value)
  message(sprintf("all - %s: %6.5f", var, pval_all))
}
# all - Population: 0.84557
# all - DNA.Source: 0.74276
# all - DNA.Extraction.Method: 0.14498
# all - MTA: 0.45966
# all - IRB: 0.60555

cat("##### Phentoype check is done\n\n")


##### (3) Duplicate checks
# check that dups are on different plates
dups$plate1 <- NA; dups$plate2 <- NA
for(i in 1:nrow(dups)) {
  dups$plate1[i] <- pm$Plate[is.element(pm$Subject_ID,dups[i,1])]
  dups$plate2[i] <- pm$Plate[is.element(pm$Subject_ID,dups[i,2])]
}
head(dups)
stopifnot(all(dups$plate1 != dups$plate2)) # T

# dups per plate
length(unique(dups$plate2)) # 35
length(unique(pm$Plate)) # 35
table(table(dups$plate2))
# 1  2 
# 6 29 
write.csv(dups, file = file.path(dir.work, "output", paste0(proj,"_dups_plate.csv")), row.names = FALSE, quote = FALSE)
cat("##### Duplicate check is done \n\n")


##### (4) Family checks 
# check that families are on same plate
pm.phen$DNA_Class <- as.numeric(pm.phen$DNA_Class)
class(pm.phen$DNA_Class)

if(!all(input.SIF$Family==0)){
  fam_check <- pm.phen %>%
    filter(DNA_Class !=3) %>%   ## added by HL: it looks like cross-dup within families is in conflict with the fule of cross-dup on different plate
    ## therefore, I got around this by limitting checks to DNA !=3
    group_by(Family) %>%
    summarise(n_plates = length(unique(Plate)))
  stopifnot(all(fam_check$n_plates == 1))
  
  ## how many families per plate?
  fams <- table(pm.phen$Family)
  fams <- fams[fams > 1]
  plt_check <- pm.phen %>%
    filter(Family %in% names(fams)) %>%
    group_by(Plate) %>%
    summarise(n_families = length(unique(Family)))
  summary(plt_check$n_families)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 11.00   28.00   30.50   31.39   35.25   41.00
} else {
  print("Not family study, family check won't be done")
}
cat("##### Family check is done\n\n")


##### (5) Cross study duplicate checks
# Cross study dups are evenly distributed across plates.
if(any(input.SIF$DNA_Class==3)){
  table(table(pm.phen$Plate[pm.phen$DNA_Class == 3]))
  # 29 plates have 1, others have 0
} else {
  print("No cross study duplicates. This check won't be done")
}
cat("##### Cross study check is done \n\n")

cat("
    ################################
    #      Bingo! You are done     #
    # Notes:",n_plates_remove,"plates were removed #
    ################################\n\n")





