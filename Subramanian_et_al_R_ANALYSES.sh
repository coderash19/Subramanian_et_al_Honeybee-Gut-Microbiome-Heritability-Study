# ====================================================================================
# CONVERTING QIIME QZA to CSV FILES FOR R
# ====================================================================================

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# README
# -----------------

## This code is to process the final .qza output files from our qiime2 pipeline and generate the OTU read count table .csv (unfiltered).
## This code does not include the filtering process, which comes after running this code and run separately on a server for faster progress
## These functions were created based on the authors file and sample naming system for this research project.

## Run step (1) once. Step (2) to (4) run 4 times with for the following 4 inputs:
#### qiime2 collapse qza files at genus level for colony 1
#### qiime2 collapse qza files at family level for colony 1
#### qiime2 collapse qza files at genus level for colony 2
#### qiime2 collapse qza files at family level for colony 2

# ------------------------------------------------------------------------------

# Step 1. Packages needed

```{r}
# Load necessary libraries
library(qiime2R)
library(tidyverse)
library(stringr)
library(purrr)
```

# ------------------------------------------------------------------------------

# Step 2. EDIT AND RUN

## ONLY FOUR THINGS TO EDIT: folder_path, the_colony, the_level, the_letter. 

```{r}

### REQUIRED INPUTS ###
#-----------------------------------------------

# Path to folder containing .qza files
folder_path <- "/path/to/Data/Microbiome/Colony2_taxa/qiime2-taxa-collapse-family-c2"

# defining blocks for the functions
the_colony = "G2"
the_level = "family"
the_letter = "f"


#-----------------------------------------------

#Inputs for the_colony:
### "G1" (this is Guts of Colony 1 bees)
### "G2" (this is Guts of Colony 2 bees)

# Inputs for the_level:
### "genus" 
### "family" 

# Inputs for the_letter:
### "g" for genus
### "f" for family
```

# ------------------------------------------------------------------------------

# 3. RUN, DO NOT EDIT - unless any other style of column ordering other than the ones below is preferred

```{r, echo=FALSE}
# Function to read .qza files and assign to global environment
read_and_assign_qza_files <- function(folder_path) {
  # List all .qza files in the specified folder
  qza_files <- list.files(path = folder_path, pattern = "\\.qza$", full.names = TRUE)
  # Function to read a .qza file and assign it to the global environment
  read_and_assign <- function(file_path) {
    # Extract file name without extension
    file_name <- str_remove(basename(file_path), "\\.qza$")
    # Read the .qza file and extract its data
    data <- read_qza(file_path)$data
    # Assign the data frame to the global environment with the file name
    assign(file_name, as.data.frame(data), envir = .GlobalEnv)
  }
  # Apply the function to all .qza files
  walk(qza_files, read_and_assign)
}


# Call the function with the specified folder path to import the qza files as data frames into global environment
read_and_assign_qza_files(folder_path)


# More definitions for functions - these things help identify and separate the files in the global environment that are needed for each of the four outputs respectively
the_pattern <- paste0("^", the_colony, ".*", the_level)
to_extract <- paste0("-collapsed-", the_level)
samples <- paste0(the_colony)
sample_NA <- paste0(the_colony, "NA")
samples_pos <- paste0(the_colony, "pos")
samples_neg <- paste0(the_colony, "neg")
final_file <- paste0(the_colony, "_", the_level)


# Combined functions to rename columns, add 'Bacteria' column, and transform data frames
process_data_frames <- function() {
  # Get the specific list of data frame names wanted from the global environment based on the data frame name
  data_frame_names <- ls(envir = .GlobalEnv, pattern = the_pattern)
  # Function to process each data frame
  process_df <- function(df_name) {
    # Get the data frame
    df <- get(df_name, envir = .GlobalEnv)
    # Extract the sample name from the data frame name
    sample_name <- str_remove(df_name, to_extract)
    # Rename columns (G1, G1NA, G1pos, G1neg)
    columns_to_rename <- c(samples, sample_NA, samples_pos, samples_neg)
    for (col_name in columns_to_rename) {
      if (col_name %in% colnames(df)) {
        colnames(df)[colnames(df) == col_name] <- sample_name
      }
    }
    # Add 'Bacteria' column containing row names
    df$Bacteria <- rownames(df)
    # Drop the original observation column (assuming it's named the same across all data frames)
    df <- df[, !names(df) %in% rownames(df)]
    # Assign the modified data frame back to the global environment
    assign(df_name, df, envir = .GlobalEnv)
  }
  # Apply the processing function to all data frames
  walk(data_frame_names, process_df)
}


# Call the combined function to process data frames
process_data_frames()


# Get the list of data frame names in the global environment
data_frame_names <- ls(pattern = the_pattern)


# Full join all data frames based on 'Bacteria' column
final_file <- data_frame_names %>%
  map(~ get(.x)) %>%
  purrr::reduce(full_join, by = "Bacteria")


# Swap the first two columns
final_file <- final_file[, c(2, 1, 3:ncol(final_file))]


# Function to standardize column names
standardize_colnames <- function(colnames_df) {
  # Identify columns to be standardized
  cols_to_standardize <- setdiff(colnames_df, "Bacteria")
  # Extract numeric suffixes and prefixes
  num_suffix <- as.numeric(gsub("^.*_(\\d+)$", "\\1", cols_to_standardize))
  prefix <- gsub(paste0("^(", the_colony, ".*)_\\d+$"), "\\1", cols_to_standardize)
  # Standardize prefixes
  standardized_prefix <- sapply(prefix, function(pref) {
    if (grepl("pos", pref, ignore.case = TRUE)) {
      return("pos")
    } else if (grepl("neg", pref, ignore.case = TRUE)) {
      return("neg")
    } else if (grepl("NA", pref, ignore.case = TRUE)) {
      return("NA")
    } else {
      return("")
    }
  })
  # Combine prefix and numeric suffix to form new column names
  new_colnames <- cols_to_standardize
  new_colnames[!grepl("Bacteria", new_colnames)] <- paste0(standardized_prefix, "_", num_suffix)
  # Combine back with "Bacteria" column
  new_colnames <- c("Bacteria", new_colnames)
  return(new_colnames)
}


# Standardize column names
colnames(final_file) <- standardize_colnames(colnames(final_file))


# Remove underscore from column names if it's at the beginning
names(final_file) <- sub("^_", "", names(final_file))


# Convert numeric column names to numeric type and reorder them
numeric_cols <- grep("^[0-9]+$", names(final_file), value = TRUE)
numeric_cols <- as.numeric(numeric_cols)  # Convert to numeric type
numeric_cols <- sort(numeric_cols)  # Sort numerically


# Combine with non-numeric column names
new_col_order <- c(setdiff(names(final_file), as.character(numeric_cols)), as.character(numeric_cols))


# Reorder columns
final_file <- final_file[new_col_order]


# Swap the first two columns
final_file <- final_file[, c(1, 40:289, 2:39)]


#manual ordering for : bactera, samples 1-250, pos1, pos2-9, pos10-18, NAs, neg1, neg2-9, neg10-18
final_file <- final_file[, c(1, 2:251, 272, 282:289, 273:281, 252:253, 254, 264:271, 255:263)]


# Rename the columns to have G1_
colnames(final_file) <- paste0(the_colony, the_letter, "_", colnames(final_file))


# Rename the first column to OTU_Bacteria"
colnames(final_file)[1] <- paste0("OTU_Bacteria")


#Check that the final file looks as expected
view(final_file)
```

# ------------------------------------------------------------------------------

# 4. SAVE FILE AS .CSV

```{r}
#Save as a csv file (not csv2)
final_name = paste0(the_colony, "_", the_level, ".csv")
write.csv(final_file, file = final_name, row.names = FALSE)
```


# ====================================================================================
# MISTAG FILTERING
# ====================================================================================

# This study was conducted with two honey bee colonies. For readability, the code below references only Colony One in folder/file names.
# Please note that this pipeline was run for Colony One samples first, and then the same code was run for Colony Two samples. 

# COLONY 1 GENUS LEVEL TAXA
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------

```{r}
/path/to/R/R-4.2.1/bin/R

library(purrr)

MaxMistagsPerBinU <- function(Samples, SpeiciesAndAbundances, MistagRate, Iters, FDR) {
  FDR <- as.numeric(as.character(FDR))
  UpperPercentile <- (1 - FDR) 
  SpeiciesAndAbundances <- SpeiciesAndAbundances[order(SpeiciesAndAbundances[,2]),]
  SeqsPerSpeciesList <- SpeiciesAndAbundances[,2]
  Taxa <- SpeiciesAndAbundances[,1]
  MistagLst <- qbinom(FDR, SeqsPerSpeciesList, MistagRate, lower.tail=F)
  MtsPerBinU <- c() 
  MinMtsPerBin <- 0
  pb <- txtProgressBar(min = 0, max = length(MistagLst), style = 3)
  for (a in 1:length(MistagLst)) {
    MaxValsU <- c()
    for (b in 1:Iters) {
      set.seed(b)
      BinLstU <- rep(0, Samples)
      du <- rdunif(MistagLst[a],Samples) }
    for (c in du) { 
      BinLstU[c] <- BinLstU[c] + 1 }
    MaxValsU[b] <- max(BinLstU, na.rm = T) 
    MtsPerBinU[a] <- quantile(MaxValsU,UpperPercentile, na.rm = T)
    if (MtsPerBinU[a] < MinMtsPerBin) {
      MtsPerBinU[a] <- MinMtsPerBin }
    MinMtsPerBin <- MtsPerBinU[a]
    setTxtProgressBar(pb, a)
  }
  df <- data.frame("Taxa"=Taxa,"SeqsPerSpecies"=SeqsPerSpeciesList,"MistagsAtFDR"=MistagLst,"MaxMTsPerBin_Uniform"=MtsPerBinU)
  rownames(df) <- df$Taxa
  return(df)
}

FilterDetections <- function(MistagThresholds, Detections){
  Out <- as.data.frame(matrix(nrow=nrow(Detections),ncol = ncol(Detections)))
  for (row in 1:nrow(Out)) {
    for (col in 1:ncol(Out)){
      Fam <- rownames(Detections)[row] 
      Out[row,col] <- ifelse(as.numeric(as.character(Detections[row,col])) <= as.numeric(as.character(MistagThresholds[Fam,4])),0,as.numeric(as.character(Detections[row,col])))
    }
  }
  rownames(Out) <- rownames(Detections) 
  colnames(Out) <- colnames(Detections)
  return(Out)
}

df <- read.csv("/path/to/R_analyses/filtering-mistags/G1_genus.csv",row.names = 1)
df[is.na(df)]=0

CMs <- colSums(df[,269:288])
Samples <- 288
TotalSeqs <- sum(colSums(df))
Abundance <- rep(rowSums(df),20)
MTs <- c()
for (i in 269:288) {
  MTs <- c(MTs,df[,i])
} 

dfLong <- data.frame("Abundance"=Abundance,"MTs"=MTs)
MistagRate <- (round(mean(CMs),0)*288)/TotalSeqs

ks.test(CMs, "pnorm", mean=mean(CMs), sd=sd(CMs)) 
CMs95CL <- mean(CMs) + qt(0.975, df = length(CMs)-1)*(sd(CMs)/length(CMs)^(1/2))
MistagRate <- (CMs95CL*Samples)/TotalSeqs 


df$Total <- rowSums(df)
SpeiciesAndAbundances <- data.frame("Taxa"=rownames(df),"Total"=df$Total)
Samples <- Samples
MistagRate <- MistagRate
Iters <- 25000
FDR <- 0.05

dfU <- MaxMistagsPerBinU(Samples,SpeiciesAndAbundances,MistagRate,Iters,FDR)
FilteredData <- FilterDetections(dfU, df[,1:268])
FilteredData = t(FilteredData) 

write.csv(FilteredData,"/path/to/R_analyses/filtering-mistags/G1_genus_filtered.csv", row.names = TRUE)
```

# COLONY 1 FAMILY LEVEL TAXA
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------

```{r}
/path/to/R/R-4.2.1/bin/R

library(purrr)

MaxMistagsPerBinU <- function(Samples, SpeiciesAndAbundances, MistagRate, Iters, FDR) {
  FDR <- as.numeric(as.character(FDR))
  UpperPercentile <- (1 - FDR) 
  SpeiciesAndAbundances <- SpeiciesAndAbundances[order(SpeiciesAndAbundances[,2]),]
  SeqsPerSpeciesList <- SpeiciesAndAbundances[,2]
  Taxa <- SpeiciesAndAbundances[,1]
  MistagLst <- qbinom(FDR, SeqsPerSpeciesList, MistagRate, lower.tail=F)
  MtsPerBinU <- c() 
  MinMtsPerBin <- 0
  pb <- txtProgressBar(min = 0, max = length(MistagLst), style = 3)
  for (a in 1:length(MistagLst)) {
    MaxValsU <- c()
    for (b in 1:Iters) {
      set.seed(b)
      BinLstU <- rep(0, Samples)
      du <- rdunif(MistagLst[a],Samples) }
    for (c in du) { 
      BinLstU[c] <- BinLstU[c] + 1 }
    MaxValsU[b] <- max(BinLstU, na.rm = T) 
    MtsPerBinU[a] <- quantile(MaxValsU,UpperPercentile, na.rm = T)
    if (MtsPerBinU[a] < MinMtsPerBin) {
      MtsPerBinU[a] <- MinMtsPerBin }
    MinMtsPerBin <- MtsPerBinU[a]
    setTxtProgressBar(pb, a)
  }
  df <- data.frame("Taxa"=Taxa,"SeqsPerSpecies"=SeqsPerSpeciesList,"MistagsAtFDR"=MistagLst,"MaxMTsPerBin_Uniform"=MtsPerBinU)
  rownames(df) <- df$Taxa
  return(df)
}

FilterDetections <- function(MistagThresholds, Detections){
  Out <- as.data.frame(matrix(nrow=nrow(Detections),ncol = ncol(Detections)))
  for (row in 1:nrow(Out)) {
    for (col in 1:ncol(Out)){
      Fam <- rownames(Detections)[row] 
      Out[row,col] <- ifelse(as.numeric(as.character(Detections[row,col])) <= as.numeric(as.character(MistagThresholds[Fam,4])),0,as.numeric(as.character(Detections[row,col])))
    }
  }
  rownames(Out) <- rownames(Detections) 
  colnames(Out) <- colnames(Detections)
  return(Out)
}

df <- read.csv("/path/to/R_analyses/filtering-mistags/G1_family.csv",row.names = 1)
df[is.na(df)]=0

CMs <- colSums(df[,269:288])
Samples <- 288
TotalSeqs <- sum(colSums(df))
Abundance <- rep(rowSums(df),20)
MTs <- c()
for (i in 269:288) {
  MTs <- c(MTs,df[,i])
} 

dfLong <- data.frame("Abundance"=Abundance,"MTs"=MTs)
MistagRate <- (round(mean(CMs),0)*288)/TotalSeqs

ks.test(CMs, "pnorm", mean=mean(CMs), sd=sd(CMs)) 
CMs95CL <- mean(CMs) + qt(0.975, df = length(CMs)-1)*(sd(CMs)/length(CMs)^(1/2))
MistagRate <- (CMs95CL*Samples)/TotalSeqs 


df$Total <- rowSums(df)
SpeiciesAndAbundances <- data.frame("Taxa"=rownames(df),"Total"=df$Total)
Samples <- Samples
MistagRate <- MistagRate
Iters <- 25000
FDR <- 0.05

dfU <- MaxMistagsPerBinU(Samples,SpeiciesAndAbundances,MistagRate,Iters,FDR)
FilteredData <- FilterDetections(dfU, df[,1:268])
FilteredData = t(FilteredData) 

write.csv(FilteredData,"/path/to/R_analyses/filtering-mistags/G1_family_filtered.csv", row.names = TRUE)
```

# ====================================================================================
# DATA FILTERING AND WRANGLING
# ====================================================================================

# STEP 1: Patriline Summary
# -------------------------------------------------------------------------------------

# Packages

```{r, echo=FALSE}
library(tidyverse)
library(vegan)
library(stringr)
```

# Load patrilines from R

```{r}
#Load Patriline info 
Colony1 = read.csv("/Users/Aishu/Desktop/YORK_UNIVERSITY/YU_School/Zayed_Research/Data/Colony1/Colony1PatrilinesforR.csv")
Colony2 = read.csv("/Users/Aishu/Desktop/YORK_UNIVERSITY/YU_School/Zayed_Research/Data/Colony2/Colony2PatrilinesforR.csv")

Summary_C1 <- Colony1 %>%
  group_by(Patriline) %>%
  summarize(No_of_offspring = n())

Summary_C2 <- Colony2 %>%
  group_by(Patriline) %>%
  summarize(No_of_offspring = n())

Summary_C1$Patriline_string = paste0("P", Summary_C1$Patriline)
Summary_C2$Patriline_string = paste0("P", Summary_C2$Patriline)
```

# Ascending order

# Colony1

```{r}

Summary_C1$Patriline_string <- factor(Summary_C1$Patriline_string, 
                                      levels = c("P0", 
                                                 Summary_C1$Patriline_string[Summary_C1$Patriline_string != "P0"][order(Summary_C1$No_of_offspring[Summary_C1$Patriline_string != "P0"])]))

# Create the plot
Colony1_plot = ggplot(Summary_C1, aes(x = Patriline_string, y = No_of_offspring)) +
  geom_bar(stat = "identity", aes(fill = ifelse(Patriline_string == "P0", "darkgrey", No_of_offspring < 5)), position = position_dodge(preserve = "single")) +  # Conditional coloring based on offspring count and adjust bar position
  scale_fill_manual(values = c("TRUE" = "lightcoral", "FALSE" = "steelblue", "darkgrey" = "darkgrey")) +  # Define colors for the condition
  labs(
    x = "Patriline",
    y = "Number of offspring"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),# Make x-axis text straight
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.ticks.y = element_line(color = "grey"),  # Add grey lines on the y-axis
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(color = "grey"),  # Set main x and y axis lines to grey
    axis.line.x.bottom = element_line(color = "grey")  # Set x-axis line color
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),  # Remove space between bars and x-axis
    breaks = scales::pretty_breaks(n = 10)  # Add increments to the y-axis
  ) +
  guides(fill = "none")  # Remove legend for fill


# Plot Colony 1
Colony1_plot

```

# Colony2

```{r}
Summary_C2$Patriline_string <- factor(Summary_C2$Patriline_string, 
                                      levels = c("P0", 
                                                 Summary_C2$Patriline_string[Summary_C2$Patriline_string != "P0"][order(Summary_C2$No_of_offspring[Summary_C2$Patriline_string != "P0"])]))


# Create the plot
Colony2_plot = ggplot(Summary_C2, aes(x = Patriline_string, y = No_of_offspring)) +
  geom_bar(stat = "identity", aes(fill = ifelse(Patriline_string == "P0", "darkgrey", No_of_offspring < 5)), position = position_dodge(preserve = "single")) +  # Conditional coloring based on offspring count and adjust bar position
  scale_fill_manual(values = c("TRUE" = "lightcoral", "FALSE" = "steelblue", "darkgrey" = "darkgrey")) +  # Define colors for the condition
  labs(
    x = "Patriline",
    y = "Number of offspring"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),# Make x-axis text straight
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.ticks.y = element_line(color = "grey"),  # Add grey lines on the y-axis
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(color = "grey"),  # Set main x and y axis lines to grey
    axis.line.x.bottom = element_line(color = "grey")  # Set x-axis line color
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),  # Remove space between bars and x-axis
    breaks = scales::pretty_breaks(n = 10)  # Add increments to the y-axis
  ) +
  guides(fill = "none")  # Remove legend for fill


# Plot Colony 2
Colony2_plot

```

# Save plots

```{r}
# Save the plot as PNG with 300 DPI
ggsave("Colony1.png", plot = Colony1_plot, width = 16, height = 6, units = "in", dpi = 300)
ggsave("Colony2.png", plot = Colony2_plot, width = 16, height = 6, units = "in", dpi = 300)


write_csv(Summary_C1, "Summary_C1.csv")
write_csv(Summary_C2, "Summary_C2.csv")

# Plots were then combined and made Publication-ready figures using Microsoft Powerpoint
```

# STEP 2: Filtering prior to Rarefaction
# -------------------------------------------------------------------------------------

## READ ME

# This is the step AFTER the mistag filtering process


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Packages

```{r, echo=FALSE}
library(tidyverse)
library(vegan)
library(stringr)
library(writexl)
library(readxl)
```

# GENUS LEVEL

## Load filtered genus level taxa collapse csvs for both colonies

```{r}
G1_genus = read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/G1_genus_silva_filtered.csv")
G2_genus = read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/G2_genus_silva_filtered.csv")
```

## Find unknown genus columns and rename those columns as unknown

```{r}

# First rename first column from X to the colony name

colnames(G1_genus)[1] <- "gut"
colnames(G2_genus)[1] <- "gut"

# Function to rename columns based on the pattern
rename_columns <- function(colnames) {
  # Apply a function to each element of 'colnames' using sapply
  sapply(seq_along(colnames), function(i) {
    # Check if the current index is 1 (first column)
    if (i == 1) {
      return(colnames[i])  # Retain the first column name and return it as is
    } else {
      colname <- colnames[i]  # Get the current column name
      # Check if the column name contains "g__"
      if (str_detect(colname, "g__")) {
        # Extract the substring starting with "g__" to the end
        genus <- str_extract(colname, "g__[^;]*$")
        # Remove "g__" from the extracted substring
        genus <- sub("g__", "", genus)
        # Check if the resulting genus name is empty
        if (genus == "") {
          return("unk_genus")  # Return "unk_genus" if the genus name is empty
        } else {
          return(genus)  # Return the extracted genus name
        }
      } else {
        return("unk_genus")  # Return "unk_genus" if "g__" is not found
      }
    }
  })
}


# G1 - genus

## Get original column names
original_colnames <- colnames(G1_genus)
## Rename columns
new_colnames <- rename_columns(original_colnames)
## Assign new column names to the data frame
colnames(G1_genus) <- new_colnames
## View the data frame to check the new column names
print(G1_genus)


# G2 - genus

## Get original column names
original_colnames <- colnames(G2_genus)
## Rename columns
new_colnames <- rename_columns(original_colnames)
## Assign new column names to the data frame
colnames(G2_genus) <- new_colnames
```

## Remove the unknown columns from each data frame

```{r}
# Function to remove columns with the name "unk_genus"
remove_unk_columns <- function(unknown) {
  # Get the column names of the data frame
  colnames <- colnames(unknown)
  
  # Logical vector to keep track of columns to retain
  retain_cols <- sapply(seq_along(colnames), function(i) {
    if (i == 1) {
      return(TRUE)  # Retain the first column
    } else {
      return(colnames[i] != "unk_genus")  # Retain columns not named "unk_genus"
    }
  })
  
  # Subset the data frame to retain only the desired columns
  unknown <- unknown[, retain_cols]
  
  return(unknown)
}


## Remove columns named "unk_genus"
G1_genus_no_unk <- remove_unk_columns(G1_genus)
G2_genus_no_unk <- remove_unk_columns(G2_genus)

```

## Remove the columns with no reads for all samples

```{r}

# Function to remove columns with all zeros
remove_zero_columns <- function(no_info) {
  # Identify columns where all values are 0
  cols_to_remove <- which(colSums(no_info == 0) == nrow(no_info))
  
  # Remove identified columns
  no_info <- no_info[, -cols_to_remove]
  
  return(no_info)
}

# G1 - genus

## Remove columns with all zeros
G1_genus_clean <- remove_zero_columns(G1_genus_no_unk)
G2_genus_clean <- remove_zero_columns(G2_genus_no_unk)

```


## Merging colony 1 and colony 2 tables based on shared columns

```{r}
#Find shared columns
shared <- intersect(names(G1_genus_clean), names(G2_genus_clean))

#Mege by shared columns
Merged_genus <- merge(G1_genus_clean, G2_genus_clean, by = shared, all = TRUE)

# Replace NA values with 0 in the merged data frame
Merged_genus[is.na(Merged_genus)] <- 0
```

## Renaming rows as sample names in prep for Vegan 

```{r}
gut_samples_values = Merged_genus$gut

rownames(Merged_genus) = gut_samples_values
```

# space

# FAMILY LEVEL

## Load filtered family level taxa collapse csvs for both colonies

```{r}
G1_family = read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/G1_family_silva_filtered.csv")
G2_family = read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/G2_family_silva_filtered.csv")
```

## Find unknown family columns and rename those columns as unknown

```{r}

# First rename first column from X to the colony name

colnames(G1_family)[1] <- "gut"
colnames(G2_family)[1] <- "gut"

# Function to rename columns based on the pattern
rename_columns <- function(colnames) {
  # Apply a function to each element of 'colnames' using sapply
  sapply(seq_along(colnames), function(i) {
    # Check if the current index is 1 (first column)
    if (i == 1) {
      return(colnames[i])  # Retain the first column name and return it as is
    } else {
      colname <- colnames[i]  # Get the current column name
      # Check if the column name contains "f__"
      if (str_detect(colname, "f__")) {
        # Extract the substring starting with "f__" to the end
        family <- str_extract(colname, "f__[^;]*$")
        # Remove "f__" from the extracted substring
        family <- sub("f__", "", family)
        # Check if the resulting family name is empty
        if (family == "") {
          return("unk_family")  # Return "unk_family" if the family name is empty
        } else {
          return(family)  # Return the extracted family name
        }
      } else {
        return("unk_family")  # Return "unk_family" if "f__" is not found
      }
    }
  })
}


# G1 - family

## Get original column names
original_colnames <- colnames(G1_family)
## Rename columns
new_colnames <- rename_columns(original_colnames)
## Assign new column names to the data frame
colnames(G1_family) <- new_colnames
## View the data frame to check the new column names
print(G1_family)

# G2 - family

## Get original column names
original_colnames <- colnames(G2_family)
## Rename columns
new_colnames <- rename_columns(original_colnames)
## Assign new column names to the data frame
colnames(G2_family) <- new_colnames
```

## Remove the unknown columns from each data frame

```{r}
# Function to remove columns with the name "unk_family"
remove_unk_columns <- function(unknown) {
  # Get the column names of the data frame
  colnames <- colnames(unknown)
  
  # Logical vector to keep track of columns to retain
  retain_cols <- sapply(seq_along(colnames), function(i) {
    if (i == 1) {
      return(TRUE)  # Retain the first column
    } else {
      return(colnames[i] != "unk_family")  # Retain columns not named "unk_family"
    }
  })
  
  # Subset the data frame to retain only the desired columns
  unknown <- unknown[, retain_cols]
  
  return(unknown)
}


## Remove columns named "unk_family"
G1_family_no_unk <- remove_unk_columns(G1_family)
G2_family_no_unk <- remove_unk_columns(G2_family)

```

## Remove the columns with no reads for all samples

```{r}

# Function to remove columns with all zeros
remove_zero_columns <- function(no_info) {
  # Identify columns where all values are 0
  cols_to_remove <- which(colSums(no_info == 0) == nrow(no_info))
  
  # Remove identified columns
  no_info <- no_info[, -cols_to_remove]
  
  return(no_info)
}


# G1 - family

## Remove columns with all zeros
G1_family_clean <- remove_zero_columns(G1_family_no_unk)
G2_family_clean <- remove_zero_columns(G2_family_no_unk)

```

## Merging colony 1 and colony 2 tables based on shared columns

```{r}
#Find shared columns
shared <- intersect(names(G1_family_clean), names(G2_family_clean))

#Merge by shared columns
Merged_family <- merge(G1_family_clean, G2_family_clean, by = shared, all = TRUE)

# Replace NA values with 0 in the merged data frame
Merged_family[is.na(Merged_family)] <- 0

```

## Renaming rows as sample names in prep for Vegan 

```{r}
gut_samples_values_2 = Merged_family$gut

rownames(Merged_family) = gut_samples_values_2
```

# Additional important data wrangling prep step for later analyses

```{r}

#Load Patriline info 
Colony1PatrilinesR = read.csv("/path/to/Data/Colony1/Colony1PatrilinesforR.csv")
Colony2PatrilinesR = read.csv("/path/to/Data/Colony2/Colony2PatrilinesforR.csv")

# ------

#Make a patrilne table 
Patrilines = as.data.frame(gut_samples_values)

Patrilines = Patrilines %>%
  dplyr:::rename(gut = gut_samples_values)

Patrilines$thorax = Patrilines$gut

#Prep patriline table to be merged with genus table
Patrilines = Patrilines %>% 
  mutate(thorax = str_replace_all(thorax, "G1g", "T1"), thorax = str_replace_all(thorax, "G2g", "T2"))

Patrilines = full_join(Patrilines, Colony1PatrilinesR, by = "thorax")
Patrilines = full_join(Patrilines, Colony2PatrilinesR, by = "thorax")

Patrilines[is.na(Patrilines)] = 0

Patrilines = Patrilines %>% mutate(Patrilines, patriline = Patriline.x + Patriline.y)

Patrilines = Patrilines %>% 
  mutate(across(c(Colony.x, Colony.y), ~str_replace_all(., "A", ""))) %>%
  mutate(Colony.x = as.numeric(Colony.x), Colony.y = as.numeric(Colony.y))

Patrilines = Patrilines %>% mutate(Patrilines, colony = Colony.x + Colony.y)

Patrilines = Patrilines %>% select(-Colony.x, -Colony.y, -Patriline.x, -Patriline.y)

Patrilines = Patrilines %>% mutate(colony = paste("A", colony, sep = ""))

Patrilines$colony[252:268] = str_replace(Patrilines$colony[252:268], "A0", "A12")
Patrilines$colony[519:536] = str_replace(Patrilines$colony[519:536], "A0", "A18")
Patrilines$colony[251] = str_replace(Patrilines$colony[251], "A0", "A12")

# ------

#create a genus + patriline info table
Patrilines_genus = full_join(Merged_genus, Patrilines, by = "gut")
Patrilines_genus = Patrilines_genus[, c(1, 95:97, 2:94)]

#Prep patriline table to be merged with family table
Patrilines = Patrilines %>% 
  mutate(gut = str_replace_all(gut, "G1g", "G1f"), gut = str_replace_all(gut, "G2g", "G2f"))

#create a family + patriline info table
Patrilines_family = full_join(Merged_family, Patrilines, by = "gut")
Patrilines_family = Patrilines_family[, c(1, 58:60, 2:57)]

#re-organise Patriline table to be more general (neither genus/family)
Patrilines = Patrilines %>% 
  mutate(gut = str_replace_all(gut, "G1f", "G1"), gut = str_replace_all(gut, "G2f", "G2"))


Patrilines_genus = Patrilines_genus %>% 
  mutate(gut = str_replace_all(gut, "G1g", "G1"), gut = str_replace_all(gut, "G2g", "G2"))

Patrilines_family = Patrilines_family %>% 
  mutate(gut = str_replace_all(gut, "G1f", "G1"), gut = str_replace_all(gut, "G2f", "G2"))

# ------

#checking percent composition of positive controls
percent_pos_genus <- Patrilines_genus %>%
  filter(grepl("pos", thorax)) %>%  # Filter rows where 'thorax' contains "pos"
  rowwise() %>%
  mutate(across(5:ncol(Patrilines_genus), ~ . / sum(c_across(5:ncol(Patrilines_genus))) * 100)) %>%
  ungroup()

percent_pos_family <- Patrilines_family %>%
  filter(grepl("pos", thorax)) %>%  # Filter rows where 'thorax' contains "pos"
  rowwise() %>%
  mutate(across(5:ncol(Patrilines_family), ~ . / sum(c_across(5:ncol(Patrilines_family))) * 100)) %>%
  ungroup()

```

## SAVE combined tables

```{r}
write_xlsx(Patrilines_genus, path = "/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_TOBErarified.xlsx")
write_xlsx(Patrilines_family, path = "/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_family_TOBErarified.xlsx")
```

# ====================================================================================
# RAREFACTION
# ====================================================================================

# STEP 1: Determine Good's Coverage 
# -------------------------------------------------------------

# Packages

```{r, echo=FALSE}
library(tidyverse)
library(vegan)
library(stringr)
library(writexl)
library(readxl)
```

```{r}
Patrilines_family_TOBErarified = read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_family_TOBErarified.csv")

Patrilines_genus_TOBErarified = read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_TOBErarified.csv")

```

# RAREFACTION PREP - determining threshold

## FAMILY
# ---------------

```{r}
# Find the number of seqs in each sample
Patrilines_family_TOBErarified$n_seqs <- rowSums(Patrilines_family_TOBErarified[, !(names(Patrilines_family_TOBErarified) %in% c("gut", "thorax", "patriline", "colony"))])
Patrilines_family_TOBErarified <- Patrilines_family_TOBErarified %>%
  select(-n_seqs, -thorax, -patriline, -colony)

# Visualizing data distribution
stats_family <- Patrilines_family_TOBErarified %>%
  rowwise() %>%
  mutate(n_seqs = sum(c_across(-gut))) %>%
  select(gut, n_seqs)

stats_family %>% 
  ggplot(aes(x=n_seqs)) +
  geom_histogram(binwidth = 1000) + 
  coord_cartesian(xlim=c(0,35000))

stats_family %>% 
  ggplot(aes(x=1, y = n_seqs)) +
  geom_jitter() + 
  scale_y_log10()

stats_family %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1:nrow(.), y=n_seqs)) +
  geom_line() + 
  coord_cartesian(xlim=c(0,100),ylim=c())

stats_family %>%
  arrange(n_seqs) %>%
  print(n=100)

# Good's coverage
stats_family_2 <- Patrilines_family_TOBErarified

stats_family_2 = pivot_longer(stats_family_2, cols = 2:ncol(stats_family_2), names_to ="family",
 values_to = "seqs")

stats_family_2 %>%
  group_by(gut) %>%
  summarize(
    n_seqs = sum(seqs), 
    n_sings = sum(seqs > 0),
    goods_family = 100*(1 - n_sings / n_seqs)) %>%
  filter(n_seqs > 0) %>%
  ggplot(aes(x=n_seqs, y=goods_family)) + 
  geom_point()

```
# CONLCUSION: use threshold of 4669 seqs, goods coverage 99.75% or higher, samples below threshold removed in .csv prior to rarefaction. New File: Patrilines_family_TOBErarified_goods.csv


## Genus
# ---------------

```{r}
# Find the number of seqs in each sample
Patrilines_genus_TOBErarified$n_seqs <- rowSums(Patrilines_genus_TOBErarified[, !(names(Patrilines_genus_TOBErarified) %in% c("gut", "thorax", "patriline", "colony"))])
Patrilines_genus_TOBErarified <- Patrilines_genus_TOBErarified %>%
  select(-n_seqs, -thorax, -patriline, -colony)

# Visualizing data distribution
stats_genus <- Patrilines_genus_TOBErarified %>%
  rowwise() %>%
  mutate(n_seqs = sum(c_across(-gut))) %>%
  select(gut, n_seqs)

stats_genus %>% 
  ggplot(aes(x=n_seqs)) +
  geom_histogram(binwidth = 1000) + 
  coord_cartesian(xlim=c())

stats_genus %>% 
  ggplot(aes(x=1, y = n_seqs)) +
  geom_jitter() + 
  scale_y_log10()

stats_genus %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1:nrow(.), y=n_seqs)) +
  geom_line() + 
  coord_cartesian(xlim=c(0,100),ylim=c())

stats_genus %>%
  arrange(n_seqs) %>%
  print(n=100)

# Good's coverage
stats_genus_2 <- Patrilines_genus_TOBErarified

stats_genus_2 = pivot_longer(stats_genus_2, cols = 2:ncol(stats_genus_2), names_to ="genus",
 values_to = "seqs")

stats_genus_2 %>%
  group_by(gut) %>%
  summarize(
    n_seqs = sum(seqs), 
    n_sings = sum(seqs > 0),
    goods_genus = 100*(1 - n_sings / n_seqs)) %>%
  filter(n_seqs > 0) %>%
  ggplot(aes(x=n_seqs, y=goods_genus)) + 
  geom_point()
```

# CONLCUSION: use threshold of 4669 seqs, goods coverage 99.75% or higher, samples below threshold removed in .csv prior to rarefaction. New File: Patrilines_genus_TOBErarified_goods.csv


# STEP 2: Rarefaction
# -------------------------------------

# 1. FAMILY

```{r}
df = na.pass(read.csv(file="/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_family_TOBErarified_goods.csv", TRUE, ",", stringsAsFactors=F, skip=0))
library(vegan)

factors = df[1:4]
observations = df[5:length(df)]

filtered = vegan::rrarefy(observations, 4669)

final = cbind(factors, filtered)

write.csv(final,"/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_family_RARIFIED.csv", row.names = FALSE)
```

# 2. GENUS

```{r}
df = na.pass(read.csv(file="/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_TOBErarified_goods.csv", TRUE, ",", stringsAsFactors=F, skip=0))
library(vegan)

factors = df[1:4]
observations = df[5:length(df)]

filtered = vegan::rrarefy(observations, 4669)

final = cbind(factors, filtered)

write.csv(final,"/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_RARIFIED.csv", row.names = FALSE)
```

# ====================================================================================
# CONVERT TO RELATIVE ABUNDANCE and PRESENCE/ABSENCE TABLES
# ====================================================================================

# STEP 1. RELATIVE ABUNDANCE
# ----------------------------------------------------------------------------------------------------------------------------------------

## GENUS

```{r}
df = na.pass(read.csv(file="/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_RARIFIED.csv", TRUE, ",", stringsAsFactors=F, skip=0))

df_factors = as.data.frame(t(df[1:4]))
df_responses = as.data.frame(t(df[5:length(df)]))

range = 1:length(df_responses)

for(i in range)
  df_responses[i] = (df_responses[i]/sum(df_responses[i]))*100

df_converted = as.data.frame(cbind(t(df_factors), t(df_responses)))

write.csv(df_converted,"/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_RARIFIED_RelAbun.csv", row.names = FALSE)
```

## FAMILY

```{r}
df2 = na.pass(read.csv(file="/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_family_RARIFIED.csv", TRUE, ",", stringsAsFactors=F, skip=0))

df2_factors = as.data.frame(t(df2[1:4]))
df2_responses = as.data.frame(t(df2[5:length(df2)]))

range = 1:length(df2_responses)

for(i in range)
  df2_responses[i] = (df2_responses[i]/sum(df2_responses[i]))*100

df2_converted = as.data.frame(cbind(t(df2_factors), t(df2_responses)))

write.csv(df2_converted,"/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_family_RARIFIED_RelAbun.csv", row.names = FALSE)
```

# STEP 2. PRESENCE/ABSENCE
# ----------------------------------------------------------------------------------------------------------------------------------------

## GENUS

```{r}
df = na.pass(read.csv(file="/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_RARIFIED.csv", TRUE, ",", stringsAsFactors=F, skip=0))

df_factors = as.data.frame(t(df[1:4]))
df_responses = as.data.frame(t(df[5:length(df)]))

range = 1:length(df_responses)

for(i in range)
  df_responses[i] = (df_responses[i]/sum(df_responses[i]))*100

df_converted = as.data.frame(cbind(t(df_factors), t(df_responses)))

write.csv(df_converted,"/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_RARIFIED_RelAbun.csv", row.names = FALSE)
```

## FAMILY

```{r}
df2 = na.pass(read.csv(file="/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_family_RARIFIED.csv", TRUE, ",", stringsAsFactors=F, skip=0))

df2_factors = as.data.frame(t(df2[1:4]))
df2_responses = as.data.frame(t(df2[5:length(df2)]))

range = 1:length(df2_responses)

for(i in range)
  df2_responses[i] = (df2_responses[i]/sum(df2_responses[i]))*100

df2_converted = as.data.frame(cbind(t(df2_factors), t(df2_responses)))

write.csv(df2_converted,"/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_family_RARIFIED_RelAbun.csv", row.names = FALSE)
```

# ====================================================================================
# GENERAL LINEAR MODELS AND PLOTS
# ====================================================================================

# Please refer to Section 2.9 Statistical Analysis in the methods section of the paper.

# For heritability of community diversity, all taxa present in a sample were used to calculate the Shannon diversity index.
# For heritability of rel-abun and pres/abs taxa that were present in at least 15% of the individuals within each colony were selected for relative abundance or presence/absence GLM analyses to ensure the data meets the assumptions of a parametric GLM model. Data was filtered down to meet this threshold on Microsoft Excel. New .csv files created were:
# Patrilines_family_RARIFIED_RelAbun_15.csv
# Patrilines_genus_RARIFIED_RelAbun_15.csv
# Patrilines_family_RARIFIED_PresAbs_15.csv
# Patrilines_genus_RARIFIED_PresAbs_15.csv


# RELATIVE ABUNDANCE
# -----------------------------------------------------------------------------------------------------------------------------------------------

# GENUS 
# -------------

```{r}
library(dplyr)
library(car)

# For patrilines with N=3

df_relabun_summary = data.frame(matrix(ncol = 10, nrow=0))
names(df_relabun_summary) = c("Bacteria", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_R2", "Col_R2", "H2")  

# Loading data
df <- read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_RARIFIED_RelAbun_15.csv")
df$patriline = as.factor(df$patriline)
df$colony = as.factor(df$colony)

# Make new column for colony_patriline and remove patrilines with n < 3
df <- df %>%
  mutate(colony_patriline = paste0(colony, "_", patriline)) %>%
  relocate(colony_patriline, .after = 4) %>%
  filter(!patriline %in% c("0", "16", "17", "18", "19", "20", "21", "22", "23"))

# remove A12_15 because its n=1 but only for A12
df <- df %>% 
  slice(-151)

df_inc = df
small_constant <- 1
df_inc[, 6:ncol(df_inc)] <- df_inc[, 6:ncol(df_inc)] + small_constant
```

### Gaussian Dist: Gilliamella, Snodgrassella
```{r}
#Gilliamella

  # 1. Run the first glm
  t1 <- glm(Gilliamella ~ colony + colony:patriline, data = df_inc, family = gaussian(link = "identity"))
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Gilliamella ~ colony, data = df_inc, family = "gaussian")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_relabun_summary = rbind(df_relabun_summary, c("Gilliamella", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))


#Snodgrassella  
  
    # 1. Run the first glm
  t1 <- glm(Snodgrassella ~ colony + colony:patriline, data = df_inc, family = "gaussian")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Snodgrassella ~ colony, data = df_inc, family = "gaussian")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_relabun_summary = rbind(df_relabun_summary, c("Snodgrassella", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))

```

### Gamma Dist: Lactobacillus, Bifidobacterium
```{r}

#Lactobacillus

  # 1. Run the first glm
  t1 <- glm(Lactobacillus ~ colony + colony:patriline, data = df_inc, family = Gamma())
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Lactobacillus ~ colony, data = df_inc, family = Gamma())
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_relabun_summary = rbind(df_relabun_summary, c("Lactobacillus", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))


#Bifidobacterium

  # 1. Run the first glm
  t1 <- glm(Bifidobacterium ~ colony + colony:patriline, data = df_inc, family = Gamma())
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Bifidobacterium ~ colony, data = df_inc, family = Gamma())
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_relabun_summary = rbind(df_relabun_summary, c("Bifidobacterium", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))

```

# FDR BH
```{r}
names(df_relabun_summary) = c("Bacteria", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  

# Convert df_relabun_summary to a data frame with the correct data types
df_relabun_summary <- as.data.frame(df_relabun_summary)
df_relabun_summary$ColPat_chi_p <- as.numeric(as.character(df_relabun_summary$ColPat_chi_p))

# Adjust pat_chi_p values using FDR (BH) correction
df_relabun_summary$ColPat_chi_p_adjust <- p.adjust(df_relabun_summary$ColPat_chi_p, method = "BH", n = length(df_relabun_summary$ColPat_chi_p))
```

#Save as CSV
```{r}
write.csv(df_relabun_summary, "/path/to/Analysis/Analyses/REVISED_Relabun_Genus_Heritability.csv", row.names = FALSE)
```

# FAMILY
# -------------

```{r}
library(dplyr)
library(car)


# For patrilines with N=3

df_relabun_summary_2 = data.frame(matrix(ncol = 10, nrow=0))
names(df_relabun_summary_2) = c("Bacteria", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_R2", "Col_R2", "H2")  

# Loading data
df <- read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_family_RARIFIED_RelAbun_15.csv")
df$patriline = as.factor(df$patriline)
df$colony = as.factor(df$colony)

# Make new column for colony_patriline
df <- df %>%
  mutate(colony_patriline = paste0(colony, "_", patriline)) %>%
  relocate(colony_patriline, .after = 4) %>%
  filter(!patriline %in% c("0", "16", "17", "18", "19", "20", "21", "22", "23"))

# remove A12_15 because its n=1
df <- df %>% 
  slice(-151)

df_inc = df
small_constant <- 1
df_inc[, 6:ncol(df_inc)] <- df_inc[, 6:ncol(df_inc)] + small_constant
```

### Gaussian Dist: Orbaceae, Neisseriaceae
```{r}
#Orbaceae

  # 1. Run the first glm
  t1 <- glm(Orbaceae ~ colony + colony:patriline, data = df_inc, family = "gaussian")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Orbaceae ~ colony, data = df_inc, family = "gaussian")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_relabun_summary_2 = rbind(df_relabun_summary_2, c("Orbaceae", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))


#Neisseriaceae  
  
    # 1. Run the first glm
  t1 <- glm(Neisseriaceae ~ colony + colony:patriline, data = df_inc, family = "gaussian")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Neisseriaceae ~ colony, data = df_inc, family = "gaussian")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_relabun_summary_2 = rbind(df_relabun_summary_2, c("Neisseriaceae", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))

```

### Gamma Dist: Lactobacillaceae, Bifidobacteriaceae
```{r}
#Lactobacillaceae

  # 1. Run the first glm
  t1 <- glm(Lactobacillaceae ~ colony + colony:patriline, data = df_inc, family = Gamma())
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Lactobacillaceae ~ colony, data = df_inc, family = Gamma())
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_relabun_summary_2 = rbind(df_relabun_summary_2, c("Lactobacillaceae", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))


#Bifidobacteriaceae

  # 1. Run the first glm
  t1 <- glm(Bifidobacteriaceae ~ colony + colony:patriline, data = df_inc, family = Gamma())
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Bifidobacteriaceae ~ colony, data = df_inc, family = Gamma())
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_relabun_summary_2 = rbind(df_relabun_summary_2, c("Bifidobacteriaceae", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))

```

# FDR BH
```{r}
names(df_relabun_summary_2) = c("Bacteria", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  

# Convert df_relabun_summary_2 to a data frame with the correct data types
df_relabun_summary_2 <- as.data.frame(df_relabun_summary_2)
df_relabun_summary_2$ColPat_chi_p <- as.numeric(as.character(df_relabun_summary_2$ColPat_chi_p))

# Adjust pat_chi_p values using FDR (BH) correction
df_relabun_summary_2$ColPat_chi_p_adjust <- p.adjust(df_relabun_summary_2$ColPat_chi_p, method = "BH", n = length(df_relabun_summary_2$ColPat_chi_p))
```

# Save as CSV
```{r}
write.csv(df_relabun_summary_2, "/path/to/Analysis/Analyses/REVISED_Relabun_Family_Heritability.csv", row.names = FALSE)
```

# PRESENCE/ABSENCE
# -----------------------------------------------------------------------------------------------------------------------------------------------

# GENUS 
# -------------
```{r}
library(dplyr)
library(car)
library(modEvA)

# For patriline with n => 3

df_presabs_summary = data.frame(matrix(ncol = 10, nrow=0))
names(df_presabs_summary) = c("Bacteria", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  

# Loading data
df <- read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_RARIFIED_PresAbs_15.csv")
df$patriline = as.factor(df$patriline)
df$colony = as.factor(df$colony)
# Make new column for colony_patriline
df <- df %>%
  mutate(colony_patriline = paste0(colony, "_", patriline)) %>%
  relocate(colony_patriline, .after = 4) %>%
  filter(!patriline %in% c("0", "16", "17", "18", "19", "20", "21", "22", "23"))

# remove A12_15 because its n=1
df <- df %>% 
  slice(-151)
```

### Binomial Dist: Commensalibacter, Snodgrassella, Frischella, Bartonella
```{r}
#Commensalibacter

  # 1. Run the first glm
  t1 <- glm(Commensalibacter ~ colony + colony:patriline, data = df, family = "binomial")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Commensalibacter ~ colony, data = df, family = "binomial")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_presabs_summary = rbind(df_presabs_summary, c("Commensalibacter", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))

#cross check with modEva
d1 = Dsquared(t1)
d2 = Dsquared(t2)

H2 = 2*(d1 - d2)
print(H2)
  

#Brachybacterium  
  
    # 1. Run the first glm
  t1 <- glm(Brachybacterium ~ colony + colony:patriline, data = df, family = "binomial")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Brachybacterium ~ colony, data = df, family = "binomial")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_presabs_summary = rbind(df_presabs_summary, c("Brachybacterium", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))


  
#Frischella  
  
    # 1. Run the first glm
  t1 <- glm(Frischella ~ colony + colony:patriline, data = df, family = "binomial")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Frischella ~ colony, data = df, family = "binomial")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_presabs_summary = rbind(df_presabs_summary, c("Frischella", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))


#Bartonella  
  
    # 1. Run the first glm
  t1 <- glm(Bartonella ~ colony + colony:patriline, data = df, family = "binomial")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Bartonella ~ colony, data = df, family = "binomial")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_presabs_summary = rbind(df_presabs_summary, c("Bartonella", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))

```

# FDR BH
```{r}
names(df_presabs_summary) = c("Bacteria", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  

# Convert df_presabs_summary to a data frame with the correct data types
df_presabs_summary <- as.data.frame(df_presabs_summary)
df_presabs_summary$ColPat_chi_p <- as.numeric(as.character(df_presabs_summary$ColPat_chi_p))

# Adjust pat_chi_p values using FDR (BH) correction
df_presabs_summary$ColPat_chi_p_adjust <- p.adjust(df_presabs_summary$ColPat_chi_p, method = "BH", n = length(df_presabs_summary$ColPat_chi_p))
```

# Save to CSV
```{r}
write.csv(df_presabs_summary, "/path/to/Analysis/Analyses/REVISED_PresAbs_Genus_Heritability.csv", row.names = FALSE)
```


# FAMILY
# -------------
 
 ```{r}
 df_presabs_summary_2 = data.frame(matrix(ncol = 10, nrow=0))
names(df_presabs_summary_2) = c("Bacteria", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  

# Loading data
df <- read.csv("/path/to//Data/Microbiome/qiime2/Post-Filter/Patrilines_family_RARIFIED_PresAbs_15.csv")
df$patriline = as.factor(df$patriline)
df$colony = as.factor(df$colony)

# Make new column for colony_patriline
df <- df %>%
  mutate(colony_patriline = paste0(colony, "_", patriline)) %>%
  relocate(colony_patriline, .after = 4) %>%
  filter(!patriline %in% c("0", "16", "17", "18", "19", "20", "21", "22", "23"))

# remove A12_15 because its n=1
df <- df %>% 
  slice(-151)
```

### Binomial Dist: Acetobacteraceae, Enterobacteriaceae, Rhizobiaceae
```{r}
#Acetobacteraceae

  # 1. Run the first glm
  t1 <- glm(Acetobacteraceae ~ colony + colony:patriline, data = df, family = "binomial")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Acetobacteraceae ~ colony, data = df, family = "binomial")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_presabs_summary_2 = rbind(df_presabs_summary_2, c("Acetobacteraceae", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))


  
#Enterobacteriaceae  
  
    # 1. Run the first glm
  t1 <- glm(Enterobacteriaceae ~ colony + colony:patriline, data = df, family = "binomial")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Enterobacteriaceae ~ colony, data = df, family = "binomial")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_presabs_summary_2 = rbind(df_presabs_summary_2, c("Enterobacteriaceae", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))

  
  
#Rhizobiaceae  
  
    # 1. Run the first glm
  t1 <- glm(Rhizobiaceae ~ colony + colony:patriline, data = df, family = "binomial")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Rhizobiaceae ~ colony, data = df, family = "binomial")
  
  # 5. Calculate deviance-based R-squared value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_presabs_summary_2 = rbind(df_presabs_summary_2, c("Rhizobiaceae", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))

```

# FDR BH
```{r}
names(df_presabs_summary_2) = c("Bacteria", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  

# Convert df_presabs_summary_2 to a data frame with the correct data types
df_presabs_summary_2 <- as.data.frame(df_presabs_summary_2)
df_presabs_summary_2$ColPat_chi_p <- as.numeric(as.character(df_presabs_summary_2$ColPat_chi_p))

# Adjust pat_chi_p values using FDR (BH) correction
df_presabs_summary_2$ColPat_chi_p_adjust <- p.adjust(df_presabs_summary_2$ColPat_chi_p, method = "BH", n = length(df_presabs_summary_2$ColPat_chi_p))
```

# Save to CSV
```{r}
write.csv(df_presabs_summary_2, "/path/to/Analysis/Analyses/Pres_Abs_Family_Heritability.csv", row.names = FALSE)
```

# ALPHA DIVERSITY
# -----------------------------------------------------------------------------------------------------------------------------------------------

```{r, echo=FALSE}
install.packages("vegan")
library(vegan)
library(viridis)
library(reshape2)
library(dplyr)
library(ggplot2)
library(car)
library(cowplot)
library(patchwork)
```

# GENUS
# ----------------

```{r}
# Load and preprocess data
dfg <- read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_RARIFIED_RelAbun.csv")
dfg$patriline <- as.character(dfg$patriline)

# Fix names for appearance
dfg <- dfg %>%
  rename("Escherichia/Shigella" = Escherichia.Shigella,
         "Candidatus Schmidhempelia" = Candidatus_Schmidhempelia,
         "Hafnia" = Hafnia.Obesumbacterium,
         "Methylorubrum" = Methylobacterium.Methylorubrum,
         "Clostridia UCG 014" = Clostridia_UCG.014)

# Filter for colonies A12 and A18 and exclude certain patrilines
dfg_A12_A18_pat <- dfg %>%
  filter(colony %in% c("A12", "A18")) %>%
  filter(!patriline %in% c("0", "16", "17", "18", "19", "20", "21", "22", "23"))

# remove A12_15 because its n=1
dfg_A12_A18_pat <- dfg_A12_A18_pat %>% 
  slice(-151)

# Calculate and check the total count of all bacteria per individual it should be 100
total_check <- dfg %>%
  rowwise() %>%
  mutate(total = sum(c_across(starts_with("Gilliamella"):ends_with("Clostridia UCG 014")), na.rm = TRUE)) %>%
  ungroup()


shannon_diversity <- diversity(dfg_A12_A18_pat[, 5:ncol(dfg_A12_A18_pat)], index = "shannon")

# Create new df with Shannon values
dfg_A12_A18_pat_Shannon = dfg_A12_A18_pat[, 1:4]
dfg_A12_A18_pat_Shannon$Shannon_G = shannon_diversity

# GLM

df_diversity_summary = data.frame(matrix(ncol = 10, nrow=0))
names(df_diversity_summary) = c("Diversity Index", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  

## Shannon Index

  # 1. Run the first glm
  t1 <- glm(Shannon_G ~ colony + colony:patriline, data = dfg_A12_A18_pat_Shannon, family = "gaussian")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2 <- Anova(t1, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi = t1_aov2[1,1]
  colony_df = t1_aov2[1,2] 
  colony_chi_p = t1_aov2[1,3]
  
  pat_chi = t1_aov2[2,1]
  pat_df = t1_aov2[2,2] 
  pat_chi_p = t1_aov2[2,3]
  
  # 4. Run second GLM
  t2 <- glm(Shannon_G ~ colony, data = dfg_A12_A18_pat_Shannon, family = "gaussian")
  
  # 5. Calculate D2 value for both GLMs: t1 and t2
  Dsquared_t1 <- 1 - (t1$deviance / t1$null.deviance)
  Dsquared_t2 <- 1 - (t2$deviance / t2$null.deviance)
  
  # 6. Calculate H2
  H2 <- 2 * (Dsquared_t1 - Dsquared_t2)
  
  # 7. Make row in table for the bacteria
  df_diversity_summary = rbind(df_diversity_summary, c("Shannon_G", colony_chi, colony_df, 
                                   colony_chi_p, pat_chi, pat_df, pat_chi_p, Dsquared_t1, Dsquared_t2, H2))
  
  names(df_diversity_summary) = c("Diversity Index", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  


# Create boxplot with patriline labels on the x-axis and no legend or colors
library(gridExtra)
  
shannong_A12 <- ggplot(subset(dfg_A12_A18_pat_Shannon, colony == "A12"), aes(x = patriline, y = Shannon_G)) +
  geom_boxplot() +
  labs(x = NULL, y = "Shannon Diversity (H')
       ") +
  coord_cartesian(ylim = c(0.5, 1.7)) +  # Set y-axis limits
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Transparent panel background
    plot.background = element_blank(),  # Transparent plot background
    #panel.grid.major = element_blank(),  # Remove all major grid lines
    #panel.grid.minor = element_blank(),  # Remove all minor grid lines
    axis.text.x = element_text(color = "black", angle = 0, size = 18),
    axis.text.y = element_text(color = "black", size = 24),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_text(size = 25, face = "bold")
  )

#Ordered
shannong_A12 <- ggplot(subset(dfg_A12_A18_pat_Shannon, colony == "A12"), 
                       aes(x = reorder_within(patriline, Shannon_G, colony, FUN = median), 
                           y = Shannon_G)) +
  geom_boxplot() +
  scale_x_reordered() +
  labs(x = "
       Patriline", y = "Shannon Diversity (H')
       ") +
  coord_cartesian(ylim = c(0.45, 1.5)) +  # Set y-axis limits
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Transparent panel background
    plot.background = element_blank(),  # Transparent plot background
    #panel.grid.major = element_blank(),  # Remove all major grid lines
    #panel.grid.minor = element_blank(),  # Remove all minor grid lines
    axis.text.x = element_text(color = "black", angle = 0, size = 18),
    axis.text.y = element_text(color = "black", size = 24),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_text(size = 25, face = "bold")
  )

shannong_A18 <- ggplot(subset(dfg_A12_A18_pat_Shannon, colony == "A18"), aes(x = patriline, y = Shannon_G)) +
  geom_boxplot() +
  labs(x = NULL, y = NULL) +  # Remove y-axis title
  coord_cartesian(ylim = c(0.5, 1.7)) +  # Set y-axis limits
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Transparent panel background
    plot.background = element_blank(),  # Transparent plot background
    #panel.grid.major = element_blank(),  # Remove all major grid lines
    #panel.grid.minor = element_blank(),  # Remove all minor grid lines
    axis.text.x = element_text(color = "black", angle = 0, size = 18),
    axis.text.y = element_text(color = "black", size = 24),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_text(size = 25, face = "bold")
  )

#Ordered
shannong_A18 <- ggplot(subset(dfg_A12_A18_pat_Shannon, colony == "A18"), 
                       aes(x = reorder_within(patriline, Shannon_G, colony, FUN = median), 
                           y = Shannon_G)) +
  geom_boxplot() +
  scale_x_reordered() +
  labs(x = "
       Patriline", y = "Shannon Diversity (H')
       ") +
  coord_cartesian(ylim = c(0.45, 1.5)) +  # Set y-axis limits
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Transparent panel background
    plot.background = element_blank(),  # Transparent plot background
    #panel.grid.major = element_blank(),  # Remove all major grid lines
    #panel.grid.minor = element_blank(),  # Remove all minor grid lines
    axis.text.x = element_text(color = "black", angle = 0, size = 18),
    axis.text.y = element_text(color = "black", size = 24),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_text(size = 25, face = "bold")
  )

# Combine the plots side by side
combinedg_plot <- grid.arrange(shannong_A12, shannong_A18, ncol = 2)


#save shannon plot
ggsave("/path/to/Analysis/Analyses/REVISED_Genus_shannon_plot.png", combinedg_plot, width = 20, height = 8, dpi = 300, bg = "transparent") 
 
```

# FAMILY
# ----------------

```{r}
# Load and preprocess data
dff <- read.csv("/path/to/Data/Microbiome/qiime2/Post-Filter/Patrilines_Family_RARIFIED_RelAbun.csv")
dff$patriline <- as.character(dff$patriline)

# Fix names for appearance
dff <- dff %>%
  rename("Clostridia UCG 014" = Clostridia_UCG.014)

# Filter for colonies A12 and A18 and exclude certain patrilines
dff_A12_A18_pat <- dff %>%
  filter(colony %in% c("A12", "A18")) %>%
  filter(!patriline %in% c("0", "16", "17", "18", "19", "20", "21", "22", "23"))

# remove A12_15 because its n=1
dff_A12_A18_pat <- dff_A12_A18_pat %>% 
  slice(-151)
  
# Calculate and check the total count of all bacteria per individual it should be 100
total_check2 <- dff %>%
  rowwise() %>%
  mutate(total = sum(c_across(starts_with("Orbaceae"):ends_with("Clostridia UCG 014")), na.rm = TRUE)) %>%
  ungroup()

shannon_diversity2 <- diversity(dff_A12_A18_pat[, 5:ncol(dff_A12_A18_pat)], index = "shannon")

# Create new df with Shannon values
dff_A12_A18_pat_Shannon = dff_A12_A18_pat[, 1:4]
dff_A12_A18_pat_Shannon$Shannon_F = shannon_diversity2


# GLM

df_diversity_summary2 = data.frame(matrix(ncol = 10, nrow=0))
names(df_diversity_summary2) = c("Diversity Index", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  


## Shannon Index

  # 1. Run the first glm
  t1.2 <- glm(Shannon_F ~ colony + colony:patriline, data = dff_A12_A18_pat_Shannon, family = "gaussian")
  
  # 2. Run the ANOVA to get the chisq value and p value
  t1_aov2.2 <- Anova(t1.2, type = 2)
  
  # 3. Define columns for results from first glm
  colony_chi.2 = t1_aov2.2[1,1]
  colony_df.2 = t1_aov2.2[1,2] 
  colony_chi_p.2 = t1_aov2.2[1,3]
  
  pat_chi.2 = t1_aov2.2[2,1]
  pat_df.2 = t1_aov2.2[2,2] 
  pat_chi_p.2 = t1_aov2.2[2,3]
  
  # 4. Run second GLM
  t2.2 <- glm(Shannon_F ~ colony, data = dff_A12_A18_pat_Shannon, family = "gaussian")
  
  # 5. Calculate R2 value for both GLMs: t1 and t2
  Dsquared_t1.2 <- 1 - (t1.2$deviance / t1.2$null.deviance)
  Dsquared_t2.2 <- 1 - (t2.2$deviance / t2.2$null.deviance)
  
  # 6. Calculate H2
  H2.2 <- 2 * (Dsquared_t1.2 - Dsquared_t2.2)
  
  # 7. Make row in table for the bacteria
  df_diversity_summary2 = rbind(df_diversity_summary2, c("Shannon_F", colony_chi.2, colony_df.2, 
                                   colony_chi_p.2, pat_chi.2, pat_df.2, pat_chi_p.2, Dsquared_t1.2, Dsquared_t2.2, H2.2))
  
  names(df_diversity_summary2) = c("Diversity Index", "Col_Chi", "Col_df", "Col_chi_p", "ColPat_Chi", "ColPat_df", "ColPat_chi_p", "ColPat_D2", "Col_D2", "H2")  


# Create boxplot with patriline labels on the x-axis and no legend or colors
library(gridExtra)

shannonf_A12 <- ggplot(subset(dff_A12_A18_pat_Shannon, colony == "A12"), aes(x = patriline, y = Shannon_F)) +
  geom_boxplot() +
  labs(x = "
       Patriline", y = "Shannon Diversity (H')
       ") +
  coord_cartesian(ylim = c(0.45, 1.5)) +  # Set y-axis limits
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Transparent panel background
    plot.background = element_blank(),  # Transparent plot background
    #panel.grid.major = element_blank(),  # Remove all major grid lines
    #panel.grid.minor = element_blank(),  # Remove all minor grid lines
    axis.text.x = element_text(color = "black", angle = 0, size = 18),
    axis.text.y = element_text(color = "black", size = 24),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_text(size = 25, face = "bold")
  )

#Ordered
shannonf_A12 <- ggplot(subset(dff_A12_A18_pat_Shannon, colony == "A12"), 
                       aes(x = reorder_within(patriline, Shannon_F, colony, FUN = median), 
                           y = Shannon_F)) +
  geom_boxplot() +
  scale_x_reordered() +
  labs(x = "
       Patriline", y = "Shannon Diversity (H')
       ") +
  coord_cartesian(ylim = c(0.45, 1.5)) +  # Set y-axis limits
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Transparent panel background
    plot.background = element_blank(),  # Transparent plot background
    #panel.grid.major = element_blank(),  # Remove all major grid lines
    #panel.grid.minor = element_blank(),  # Remove all minor grid lines
    axis.text.x = element_text(color = "black", angle = 0, size = 18),
    axis.text.y = element_text(color = "black", size = 24),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_text(size = 25, face = "bold")
  )

print(shannonf_A12)

shannonf_A18 <- ggplot(subset(dff_A12_A18_pat_Shannon, colony == "A18"), aes(x = patriline, y = Shannon_F)) +
  geom_boxplot() +
  labs(x = "
       Patriline", y = NULL) +  # Remove y-axis title
  coord_cartesian(ylim = c(0.45, 1.5)) +  # Set y-axis limits
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Transparent panel background
    plot.background = element_blank(),  # Transparent plot background
    #panel.grid.major = element_blank(),  # Remove all major grid lines
    #panel.grid.minor = element_blank(),  # Remove all minor grid lines
    axis.text.x = element_text(color = "black", angle = 0, size = 18),
    axis.text.y = element_text(color = "black", size = 24),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_text(size = 25, face = "bold")
  )

print(shannonf_A18)

#Ordered
shannonf_A18 <- ggplot(subset(dff_A12_A18_pat_Shannon, colony == "A18"), 
                       aes(x = reorder_within(patriline, Shannon_F, colony, FUN = median), 
                           y = Shannon_F)) +
  geom_boxplot() +
  scale_x_reordered() + 
  labs(x = "
       Patriline", y = "Shannon Diversity (H')
       ") +
  coord_cartesian(ylim = c(0.45, 1.5)) +  # Set y-axis limits
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Transparent panel background
    plot.background = element_blank(),  # Transparent plot background
    #panel.grid.major = element_blank(),  # Remove all major grid lines
    #panel.grid.minor = element_blank(),  # Remove all minor grid lines
    axis.text.x = element_text(color = "black", angle = 0, size = 18),
    axis.text.y = element_text(color = "black", size = 24),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_text(size = 25, face = "bold")
  )

# Combine the plots side by side
combinedf_plot <- grid.arrange(shannonf_A12, shannonf_A18, ncol = 2)


#save shannon plot
ggsave("/path/to/Analysis/Analyses/REVISED_Family_shannon_plot.png", combinedf_plot, width = 20, height = 8, dpi = 300, bg = "transparent") 
```

# BETA DIVERSITY
# ------------------------------------------------------------------------------------------------------------------------------------------------

```{r}

# PACKAGES

library(vegan)
library(ape)
library(viridis)
library(reshape2)
library(dplyr)
library(ggplot2)
library(car)
library(cowplot)
library(patchwork)
library(tidyr)
library(tidyverse)
```

# PREP WORK FOR PCOA

## GENUS

```{r}
# Load and preprocess data
dfg <- read.csv("/Users/Aishu/Desktop/YORK_UNIVERSITY/YU_School/Zayed_Research/Data/Microbiome/qiime2/Post-Filter/Patrilines_genus_RARIFIED_RelAbun.csv")

dfg$patriline <- as.character(dfg$patriline)

# Fix names for appearance
dfg <- dfg %>%
  rename("Escherichia/Shigella" = Escherichia.Shigella,
         "Candidatus Schmidhempelia" = Candidatus_Schmidhempelia,
         "Hafnia" = Hafnia.Obesumbacterium,
         "Methylorubrum" = Methylobacterium.Methylorubrum,
         "Clostridia UCG 014" = Clostridia_UCG.014)

# Filter for colonies A12 and A18 and exclude certain patrilines
dfg_A12_A18_pat <- dfg %>%
  filter(colony %in% c("A12", "A18")) %>%
  filter(!patriline %in% c("0", "16", "17", "18", "19", "20", "21", "22", "23"))

# Remove A12_15 because its n=1
dfg_A12_A18_pat <- dfg_A12_A18_pat %>% 
  slice(-151)

# Add 'A' and 'B' in front of Patrilines
dfg_A12_A18_pat$patriline <- factor(
  paste0(
    ifelse(dfg_A12_A18_pat$colony == "A12", "A", "B"), 
    dfg_A12_A18_pat$patriline
  ), 
  levels = c(paste0("A", 1:14), paste0("B", 1:15))
)


# Calculate Bray-Curtis dissimilarity
bray_curtis <- vegdist(dfg_A12_A18_pat[, 5:ncol(dfg_A12_A18_pat)], method = "bray")

# Perform PCoA
pcoa_bray <- cmdscale(bray_curtis, eig = TRUE, k = 2)

# Extract PCoA coordinates
pcoa_coords <- as.data.frame(pcoa_bray$points)
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")

# Add colony and patriline information
pcoa_coords$colony <- dfg_A12_A18_pat$colony
pcoa_coords$patriline <- dfg_A12_A18_pat$patriline
```

## FAMILY

```{r}
# Load and preprocess data
dff <- read.csv("/Users/Aishu/Desktop/YORK_UNIVERSITY/YU_School/Zayed_Research/Data/Microbiome/qiime2/Post-Filter/Patrilines_Family_RARIFIED_RelAbun.csv")

dff$patriline <- as.character(dff$patriline)

# Fix names for appearance
dff <- dff %>%
  rename("Clostridia UCG 014" = Clostridia_UCG.014)

# Filter for colonies A12 and A18 and exclude certain patrilines
dff_A12_A18_pat <- dff %>%
  filter(colony %in% c("A12", "A18")) %>%
  filter(!patriline %in% c("0", "16", "17", "18", "19", "20", "21", "22", "23"))

# Remove A12_15 because its n=1
dff_A12_A18_pat <- dff_A12_A18_pat %>% 
  slice(-151)

# Add 'A' and 'B' in front of Patrilines
dff_A12_A18_pat$patriline <- factor(
  paste0(
    ifelse(dff_A12_A18_pat$colony == "A12", "A", "B"), 
    dff_A12_A18_pat$patriline
  ), 
  levels = c(paste0("A", 1:14), paste0("B", 1:15))
)


# Calculate Bray-Curtis dissimilarity
bray_curtis2 <- vegdist(dff_A12_A18_pat[, 5:ncol(dff_A12_A18_pat)], method = "bray")

# Perform PCoA
pcoa_bray2 <- cmdscale(bray_curtis2, eig = TRUE, k = 2)

# Extract PCoA coordinates
pcoa_coords2 <- as.data.frame(pcoa_bray2$points)
colnames(pcoa_coords2) <- c("PCoA1", "PCoA2")

# Add colony and patriline information
pcoa_coords2$colony <- dff_A12_A18_pat$colony
pcoa_coords2$patriline <- dff_A12_A18_pat$patriline

```

# GETTING AXIS FOR PCOAs

```{r}
# Revised method for calculating variance explained
summarize_axes <- function(pcoa_object, dataset) {
  eig <- pcoa_object$eig
  n_axes <- min(10, length(eig))
  pct <- 100 * eig[seq_len(n_axes)] / sum(eig[eig > 0])
  
  data.frame(
    Dataset = dataset,
    Axis = paste0("PCoA", seq_len(n_axes)),
    Eigenvalue = eig[seq_len(n_axes)],
    Percent_explained = pct,
    Cumulative_percent = cumsum(pct)
  )
}

summarize_axes <- function(pcoa_object, dataset) {
  eig <- pcoa_object$values$Eigenvalues
  n_axes <- min(10, length(eig))
  pct <- 100 * eig[seq_len(n_axes)] / sum(eig[eig > 0])

  data.frame(Dataset = dataset, Axis = paste0("PCoA", seq_len(n_axes)), Eigenvalue = eig[seq_len(n_axes)], Percent_explained = pct, Cumulative_percent = cumsum(pct))
}


# Get eigenvalues using ape::pcoa
pcoa_bray_values <- pcoa(bray_curtis)
pcoa_bray2_values <- pcoa(bray_curtis2)


# Calculate variance explained
variance_table <- rbind(
  summarize_axes(pcoa_bray_values, "Genus"),
  summarize_axes(pcoa_bray2_values, "Family")
)

print(variance_table)


# Get PC1 and PC2 values for Genus
pc_g <- variance_table %>%
  filter(Dataset == "Genus", Axis %in% c("PCoA1", "PCoA2")) %>%
  pull(Percent_explained)


# Get PC1 and PC2 values for Family
pc_f <- variance_table %>%
  filter(Dataset == "Family", Axis %in% c("PCoA1", "PCoA2")) %>%
  pull(Percent_explained)


# Print values
pc_g
pc_f

# Create axis labels
genus_x <- paste0("PC1 (", round(pc_g[1], 2), "% variance explained)")
genus_y <- paste0("PC2 (", round(pc_g[2], 2), "% variance explained)")

family_x <- paste0("PC1 (", round(pc_f[1], 2), "% variance explained)")
family_y <- paste0("PC2 (", round(pc_f[2], 2), "% variance explained)")

```

# PLOTTING PCOAS
```{r}
# COLONY 1

# Balanced mix of open and filled shapes
custom_shapes <- c(16, 15, 17, 0, 18, 1, 19, 2, 20, 5, 21, 6, 22, 23)

# Colors using the "viridis" palette
custom_colors <- viridis_pal(option = "turbo")(14)


# Plot for colony A12 - GENUS
plotg_A12 <- ggplot(pcoa_coords %>% filter(colony == "A12"), aes(x = PCoA1, y = PCoA2, color = patriline, shape = patriline)) +
  geom_point(size = 8) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  labs(x = genus_x,
       y = genus_y) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_text(size = 50, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 30),
    legend.position = "none",
    legend.title = element_text(size = 40, face = "bold"),
    legend.text = element_text(size = 35),
    legend.key.size = unit(2, "cm")
  )


# Plot for colony A12 - FAMILY
plotf_A12 <- ggplot(pcoa_coords2 %>% filter(colony == "A12"), aes(x = PCoA1, y = PCoA2, color = patriline, shape = patriline)) +
  geom_point(size = 8) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  labs(x = family_x,
       y = family_y) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_text(size = 50, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 30),
    legend.position = "right",
    legend.title = element_text(size = 40, face = "bold"),
    legend.text = element_text(size = 35),
    legend.key.size = unit(2, "cm")
  )


# Combine A12 plots
PCoA_A12_plot <- (plotg_A12 + plotf_A12 + plot_layout(guides = "collect")) &
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

PCoA_A12_plot


# COLONY 2

# Balanced mix of open and filled shapes
custom_shapes2 <- c(15, 16, 0, 17, 1, 18, 2, 19, 5, 20, 6, 21, 22, 23, 24)

# Colors using the "magma" palette
custom_colors2 <- viridis_pal(option = "magma", begin = 0, end = 0.75)(15)


# Plot for colony A18 - GENUS
plotg_A18 <- ggplot(pcoa_coords %>% filter(colony == "A18"), aes(x = PCoA1, y = PCoA2, color = patriline, shape = patriline)) +
  geom_point(size = 8) +
  scale_shape_manual(values = custom_shapes2) +
  scale_color_manual(values = custom_colors2) +
  labs(x = genus_x,
       y = genus_y) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_text(size = 50, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 30),
    legend.position = "right",
    legend.title = element_text(size = 40, face = "bold"),
    legend.text = element_text(size = 35),
    legend.key.size = unit(2, "cm")
  )


# Plot for colony A18 - FAMILY
plotf_A18 <- ggplot(pcoa_coords2 %>% filter(colony == "A18"), aes(x = PCoA1, y = PCoA2, color = patriline, shape = patriline)) +
  geom_point(size = 8) +
  scale_shape_manual(values = custom_shapes2) +
  scale_color_manual(values = custom_colors2) +
  labs(x = family_x,
       y = family_y) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_text(size = 50, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 30),
    legend.position = "right",
    legend.title = element_text(size = 40, face = "bold"),
    legend.text = element_text(size = 35),
    legend.key.size = unit(2, "cm")
  )


# Combine A18 plots
PCoA_A18_plot <- (plotg_A18 + plotf_A18 + plot_layout(guides = "collect")) &
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

PCoA_A18_plot

```

# SAVE REVISED PLOTS

```{r}
ggsave("REVISED-MANUSCRIPT_PCoA_A12.tiff", plot = PCoA_A12_plot, width = 33, height = 13, units = "in", dpi = 600, compression = "lzw")

ggsave("REVISED-MANUSCRIPT_PCoA_A18.tiff", plot = PCoA_A18_plot, width = 33, height = 13, units = "in", dpi = 600, compression = "lzw")
```


