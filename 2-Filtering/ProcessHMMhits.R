# Load required libraries
pacman::p_load(tidyverse)

# Load the clean HMM score data
dfFull <- read_delim("HMMERhits.tblout", delim = "\t", col_names = c("Identifier", "Family","HMMscore"))

# Make a new data frame with data needed
df <- tibble(Identifier=unique(dfFull$Identifier))
families <- c(unique(dfFull$Family))
df[,families] <- NA

# Add HMM score from dfFull to dfScores
dfScores <- df
for(i in 1:nrow(dfFull)) {
  familyColumn <- which(names(dfScores) == dfFull$Family[i])
  idRow <- which(dfScores$Identifier == dfFull$Identifier[i])
  dfScores[idRow,familyColumn] <- dfFull$HMMscore[i]
}

# Get the top two HMM score values and their families 
dfScores2 <- dfScores
dfScores2[,c("FirstValue","SecondValue")] <- NA
for(i in 1:nrow(dfScores2)) {
  x <- na.omit(as.numeric(dfScores[i,2:15]))
  dfScores2$FirstValue[i] <- max(x)
  dfScores2$SecondValue[i] <- max(x[x != max(x)])
}
dfScores2$FirstHit <- sapply(1:nrow(dfScores2), function(x) colnames(dfScores2)[which(dfScores2[x,2:15] == dfScores2$FirstValue[x])+1])
dfScores2$SecondHit <- sapply(1:nrow(dfScores2), function(x) colnames(dfScores2)[which(dfScores2[x,2:15] == dfScores2$SecondValue[x])+1])

# Calculate 'fold change' between first and second hits
dfScores3 <- dfScores2
dfScores3$SecondValue[dfScores3$SecondValue == -Inf] <- 1.0
dfScores3$SecondHit[dfScores3$SecondHit == "character(0)"] <- "None"
dfScores3$FoldDifference <- dfScores3$FirstValue/dfScores3$SecondValue
dfScores3$AbsoluteDifference <- dfScores3$FirstValue-dfScores3$SecondValue
dataProcessed <- dfScores3 %>% rowwise() %>% mutate(FirstHit = paste(FirstHit, collapse=',')) %>% 
  mutate(SecondHit = paste(SecondHit, collapse=',')) %>% ungroup()

# Filtering and finalize hits
dataFinal <- dataProcessed %>% mutate(FinalHit = case_when(
  FirstValue > median(dataProcessed$FirstValue) & SecondValue == 1  ~ FirstHit,
  FirstValue <= median(dataProcessed$FirstValue)  ~ "None",
  FoldDifference > 1.5 & AbsoluteDifference > median(dataProcessed$FirstValue)  ~ FirstHit,
  .default = "Multiple" 
  ))

# Remove unwanted 'FinalHit' categories and split identifiers/accessions
dataFinal2 <- dataFinal %>% filter(FinalHit != "AddA,AdnB" & FinalHit != "Multiple" & FinalHit != "None")
ggplot(dataFinal2, aes(x=FirstValue)) + geom_histogram(color="darkblue", fill="lightblue") + facet_wrap(~ FinalHit, ncol=5, scales = 'free')
dataFinal2 <- dataFinal2 %>% select(Identifier, FinalHit)
dataFinal2 <- dataFinal2 %>% separate_wider_delim(Identifier, delim = "-", names = c("PhyClaOrd", "Genus", "AssemblyAccession", "ProteinID"), cols_remove = F)

# Write final data to file
write.table(dataFinal2, file = "ProcessedHMMhits.tsv", quote = F, sep = "\t", row.names = F, col.names = T)


