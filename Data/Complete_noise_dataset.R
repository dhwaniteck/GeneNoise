rna_noise_acar <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_acar.csv", header = TRUE, sep = ",", dec = ".")
rna_noise_gasch <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_gasch.csv", header = TRUE, sep = ",", dec = ".")
rna_noise_huang <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_huang.csv", header = TRUE, sep = ",", dec = ".")
rna_noise_jariani <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_jariani.csv", header = TRUE, sep = ",", dec = ".")
rna_noise_steinmetz <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_steinmetz.csv", header = TRUE, sep = ",", dec = ".")
rna_noise_steinmetz$X <- NULL

rna_noise <- merge(rna_noise_acar,rna_noise_gasch,by="ORF",all.x = TRUE, all.y=TRUE)
rna_noise <- merge(rna_noise,rna_noise_huang,by="ORF",all.x = TRUE, all.y=TRUE)
rna_noise <- merge(rna_noise,rna_noise_jariani,by="ORF",all.x = TRUE, all.y=TRUE)
rna_noise <- merge(rna_noise,rna_noise_steinmetz,by="ORF",all.x = TRUE, all.y=TRUE)

promoter_sequences <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/yeast_promoters.csv", header = TRUE, sep = ",", dec = ".")

rna_noise <- merge(rna_noise,promoter_sequences,by="ORF")

colnames(rna_noise)

write.csv(rna_noise,"/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_complete.csv")

# Load necessary libraries
library(ggplot2)
library(patchwork) # Use patchwork for grid layout

# Assume your data is in a data frame called "data"

# Extract unique conditions
conditions <- c("DMSO", "gua", "MPA", "guaMPA", "stressed", "unstressed", 
                "2hr", "16hr", "36hr", "glu6hr", "glu12hr", 
                "lag1hr", "lag3hr", "glumal", "YPD", "YJM")

# Initialize a list to store plots
plot_list <- list()

# Loop through conditions and create plots
for (condition in conditions) {
  # Dynamically construct column names
  dm_col <- paste0("DM_", condition)
  avgexp_col <- paste0("AVGEXP_", condition)
  
  # Check if columns exist in the data
  if (dm_col %in% names(rna_noise) & avgexp_col %in% names(rna_noise)) {
    # Create the scatter plot
    p <- ggplot(rna_noise, aes_string(x = avgexp_col, y = dm_col)) +
      geom_point() +
      scale_x_log10() + # Apply log transformation to the x-axis
      labs(title = condition, x = paste("AVGEXP_", condition, sep = ""), y = paste("DM_", condition, sep = "")) +
      theme_minimal()
    
    # Add plot to the list
    plot_list[[condition]] <- p
  } else {
    cat("Skipping", condition, "as columns are missing.\n")
  }
}

# Combine all plots in a grid
all_plots <- wrap_plots(plot_list, ncol = 4) # Adjust ncol as needed
all_plots

plots <- list()

# Loop through all conditions to create individual plots
for (condition in conditions) {
  avgexp_col <- paste0("AVGEXP_", condition)
  cv_col <- paste0("CV_", condition)
  
  # Generate the plot for the current condition
  p <- ggplot(rna_noise, aes(
    x = log(!!sym(avgexp_col)), 
    y = log((!!sym(cv_col))^2)
  )) +
    geom_point(alpha = 0.6) +
    labs(
      title = condition,
      x = "Log(AVGEXP)",
      y = "Log(CV^2)"
    ) +
    theme_minimal()
  
  # Append plot to the list
  plots[[condition]] <- p
}

# Combine all plots into a grid
grid <- wrap_plots(plots, ncol = 4) # Adjust ncol for desired layout
print(grid)

jariani <- rna_noise[,c("ORF","DM_glu6hr","promoter")]
jariani <- jariani[!is.na(jariani$DM_glu6hr) & !is.na(jariani$promoter),]
jariani <- jariani[-which(nchar(jariani$promoter)<1000),]
table(nchar(jariani$promoter))

hist(jariani$DM_glu6hr)
table(jariani$DM_glu6hr>0.01)

write.csv(jariani,"/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/glu6hr_CNNdata.csv")

# Doing Fisher's exact test on noisy genes from all conditions
library(dplyr)

# Example data frame
# Replace this with your actual data
# noise <- read.csv("your_data.csv")

# Step 1: Select columns with "DM_" in their names
dm_columns <- rna_noise %>% select(starts_with("DM_"))

# Step 3: Perform Fisher's Exact Test for each pair of DM_ columns
results <- list()

for (col1 in names(dm_columns)) {
  for (col2 in names(dm_columns)) {
    if (col1 != col2) {  # Avoid testing a column against itself
      # Create a contingency table
      contingency_table <- table(dm_columns[[col1]] > 0.01, dm_columns[[col2]] > 0.01)
      
      # Check if the table is 2x2
      if (all(dim(contingency_table) == c(2, 2))) {
        # Perform Fisher's Exact Test
        fisher_result <- fisher.test(contingency_table)
        
        # Save results
        results[[length(results) + 1]] <- data.frame(
          Column1 = col1,
          Column2 = col2,
          OddsRatio = fisher_result$estimate,
          PValue = fisher_result$p.value
        )
      }
    }
  }
}

# Step 4: Combine results into a single data frame
results_df <- do.call(rbind, results)

# View results
results_df <- results_df %>% arrange(PValue)
print(results_df)
