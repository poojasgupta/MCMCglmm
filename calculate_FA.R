# Load necessary libraries
library(lme4)
library(dplyr)
library(tidyr)

# Load your dataset (modify the path as needed)
data <- read.csv("example_FA_data.csv")

# Select relevant columns for left and right tarsus measurements, species, and Major.Loc (trapping site)
# Similar approach can be used for wing measurements
rData <- data %>%
  select(SN, Species, Major.Loc.x, Left.T1, Left.T2, Left.T3, Right.T1, Right.T2, Right.T3)

# Reshape data into long format for left-side tarsus measurements
left_data <- rData %>%
  pivot_longer(cols = c(Left.T1, Left.T2, Left.T3), 
               names_to = 'Trait', values_to = 'Measure') %>%
  mutate(Side = 1)  # 1 represents left side

# Reshape data into long format for right-side tarsus measurements
right_data <- rData %>%
  pivot_longer(cols = c(Right.T1, Right.T2, Right.T3), 
               names_to = 'Trait', values_to = 'Measure') %>%
  mutate(Side = -1)  # -1 represents right side

# Combine left and right tarsus data into one dataset
FA_Data <- bind_rows(left_data, right_data)

# Remove redundant columns (Left.T1, Right.T1, etc.)
FA_Data <- FA_Data %>%
  select(SN, Species, Major.Loc.x, Trait, Measure, Side)

# Sort the data by individual and side for consistency
FA_Data <- FA_Data %>%
  arrange(SN, Side, Trait)

# Fit a linear mixed-effects model to calculate FA for tarsus, including Species and Major.Loc as fixed effects
mod <- lmer(Measure ~ Side + Species + Major.Loc.x + (Side | SN), data = FA_Data)

# Extract random effects (individual asymmetry effects)
RandEff <- ranef(mod)$SN

# Create a dataframe with calculated FA values for each individual
FA_results <- data.frame(SN = as.numeric(row.names(RandEff)),
                         FA.T = abs(RandEff[, 2]))  # Absolute value of the asymmetry effect

# Save the results to a CSV file
# write.csv(FA_results, "FA_results_tarsus.csv", row.names = FALSE)

# Merge FA results back into the original data
merged_data <- left_join(data, FA_results, by = "SN")

# Save the merged dataset to a CSV file
write.csv(merged_data, "merged_data_with_FA.csv", row.names = FALSE)
