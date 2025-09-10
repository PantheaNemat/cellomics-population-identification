## Analysis single cell data derived from Columbus

# determine cut-offs based on PCA analyses and different clustering algorithms of all your single cell values
# instead of manual inspection in Columbus
# let the data speak for itself :) 

# I have used three different variables to identify neuronal nuclei 
# 1) Hoechst intensity (too little Hoechst intensity -> not a nucleus, too high Hoechst intensity -> more likely an apoptopic cell)
# 2) Roundness of identified nuclei (nucleus needs to be round, otherwise might be a different structure or dividing cell)
# 3) Size (too small -> most likely apoptopic cell, too large -> different cell type, dividing cell etc.)
# optional: 4) MAP2 intensity - in case of primary cultures (nuclei should be MAP2 positive)

# Identified nuclei are then further divided into GFP positive or negative nuclei (or even more subdivisions)
# based on GFP intensity using hierarchical clustering

# data is then ready to be analysed using a statistical test of your choice 
# given the dependency of the observations within and across well plates
# I would recommend multilevel analyses which is also included in this script 


# let's start!

# download and load packages by running packages script

# Clear existing work space objects 
rm(list = ls())
graphics.off()

#create a folder with all the files containing single cell data 

# set a working directory to this folder
setwd("~/Documents/Drebrin A/SINEUP/Export N2A Cells/Dbn1_N2A_SINEUP_sc_only")

# create a list with all the file names 
file_list <- list.files(path ="~/Documents/Drebrin A/SINEUP/Export N2A Cells/Dbn1_N2A_SINEUP_sc_only")
# initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable

rawdata <- data.frame()
df <- data.frame()

# for loop will loop through folder, read in data file and attach it to existing data frame 
# result is a single big data file containing ALL your daya
for (i in 1:length(file_list)) {
  temp_data <- read_csv(file_list[i], col_names = TRUE) 
  rawdata <- rbind(rawdata, temp_data) 
}


# check data integrity
# check whether all your columns are in, all your files etc
head(rawdata)
tail(rawdata)
dim(rawdata)
names(rawdata)
summary(rawdata)

df <- rawdata

# set working directory to folder where you want to store big data file 
setwd("~/Documents/Drebrin A/SINEUP/N2A Cells/scDATA")

# store big data file 
write.table(rawdata,
            file = "SINEUP_Dbn1_N2ACells_scDATA_original.csv",
            sep = ";", dec = ",", row.names = F
)

# I would recommend renaming the variable names to something shorter/more handy 
# variable names in columbus are usually quite long and entail spaces which can be annoying 

#### #### ####
# once you have big data file, you can simply run it from here 
df <- read.table("SINEUP_Dbn1_N2ACells_scDATA.csv", sep=",", dec=".", header=TRUE)

# if you want to perform multilevel analysis later and you scanned several well plates
# you can create a unique identifier by pasting the well plate name (screen name) and the well name 
df$WellPlate_Well <- paste(df$ScreenName, df$WellName, sep = "_")


## create variable containing information on condition 
# in case you still need it 
df$Condition <- 0

# adapt to your own well plate set up 
df <- df %>%
  mutate(Condition = case_when(
    Column == 8 ~ "CTL",
    Column == 5 ~ "SU1",
    Column == 6 ~ "SU2",
    Column == 7 ~ "SU3",
    
    TRUE ~ NA_character_
  ))

# check whether naming worked out correctly 
combinations <- unique(df[, c("ScreenName", "WellName", "Condition")]) 


## get overview of your data 
hist(df$Dbn1_MeanInten, breaks = 500)
hist(df$GFP_MeanInten, breaks = 500)
hist(df$NuclRoundness, breaks = 500)
hist(df$NucleusArea, breaks = 500)
hist(df$Hoechst_MeanInten, breaks = 500)

psych::describe(df)


## turn all your variables that are factors into true factors 
# also only more important for later multilevel analysis 
df$ScreenName <- as.factor(df$ScreenName)
df$WellName <- as.factor(df$WellName)
df$WellPlate_Well <- as.factor(df$WellPlate_Well)
df$Condition <- as.factor(df$Condition)


## data file should now be ready to determine cut offs for your variables

#####
# to get a feeling for your data and which variables seem to explain variation in your data
# it might be helpful to first run a PCA

# select all variables that might aid in identifying neurons into a data frame
# I would not put the staining intensity of your protein of interest here (in my case Drebrin), it should really be about the morphology/neuron identification
df_pca <- df[, c(7:9)] 
 
# Identify principle components and loading of each parameter on each component
res.pca <- prcomp(df_pca, scale = TRUE)
print(res.pca)

# Identify and plot how much variance is explained by each dimension 
eig.val <- get_eigenvalue(res.pca)
eig.val

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

var <- get_pca_var(res.pca)
var

# Coordinates
head(var$coord)
head(var$cor)
# Cos2: quality on the factor map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

# Coordinates of variables
head(var$coord, 4)

# Visualisation of parameters on each dimension
fviz_pca_var(res.pca, col.var = "black")
corrplot(var$cos2, is.corr = FALSE)

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(res.pca,
             col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# Contributions of variables to PC1+2
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)


# Show clustering according to number of identified clusters
fviz_pca_biplot(res.pca,
                geom.ind = "point", # show points only (nbut not "text")
                addEllipses = FALSE, # Concentration ellipses
                legend.title = "Cluster",
                pointsize = 2
)

# Show clustering according to number of identified clusters
fviz_pca_biplot(res.pca,
                geom.ind = "point", # show points only (nbut not "text")
                col.ind = df$Condition, # color by cluster
                palette = c("#FF0000", "#FFFF00", "#00FF00", "#0000FF"), # Add as many color codes as determined clusters
                addEllipses = FALSE, # Concentration ellipses
                legend.title = "Cluster",
                pointsize = 2
)


# arrows that are perpendicular to each indicate that variables are not correlated with each other 
# the closer the arrows the more positively correlated the two variables
# arrows in opposing directions indicate variables are negatively correlated with each other



#### Identify neuronal nuclei using Hoechst intensity, roundness and area
## Determine cut-offs for each variable by performing PCA                              

## 1.) Hoechst intensity 

## Hierarchical clustering
# Identification of clusters and visualisation 
# select column containing Hoechst intensity 
df_pca <- df[, c("Hoechst_MeanInten")]

df_pca <- as.data.frame(df_pca)


## K-means clustering
set.seed(123) # always set specific seed in order to be able to replicate results

# Elbow approach
# Function to compute total within-cluster sum of square
wss <- function(k) {
  kmeans(df_pca, k, nstart = 10)$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15
# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares"
)
# determine perfect amount of clusters based on most extreme bend in plot 

# specify in centers how many clusters you have identified in the plot using the Elbow approach
km <- kmeans(df$Hoechst_MeanInten,centers=2)
df$clust_hoechst <- as.factor(km$cluster)

# plot your data using the determined classification
ggplot(df, aes(x=Hoechst_MeanInten)) + 
  geom_histogram(aes(fill=clust_hoechst),
                 binwidth=0.5)+
  stat_density(geom="line", color="red")

# in example, no clear segregation of two distributions visible, can also choose not to use this variable if you are not convinced 
# also play around with bigger/smaller binwidth/bins, sometimes then you start to see multiple distributions 


# 2.) Area
## Hierarchical clustering
# Identification and visualisation 
df_pca <- df[, c("NucleusArea")]

df_pca <- as.data.frame(df_pca)

## K-means clustering
set.seed(123) # always set specific seed in order to be able to replicate results

# Elbow approach
# Function to compute total within-cluster sum of square
wss <- function(k) {
  kmeans(df_pca, k, nstart = 10)$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15
# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares"
)
# determine perfect amount of clusters based on most extreme bend in plot 

# specify in centers how many clusters you have identified in the plot using the Elbow approach
km <- kmeans(df$NucleusArea,centers=2)
df$clust_area <- as.factor(km$cluster)

# plot your data using the determined classification
ggplot(df, aes(x=NucleusArea)) + 
  geom_histogram(aes(fill=clust_area),
                 binwidth=0.5)+
  stat_density(geom="line", color="red")

# in example two distributions visible - normal distribution + right flank 
# larger nuclei most nuclei that are dividing, so choose cluster of smaller nuclei here 



## 3.) Roundness
## Hierarchical clustering
# Identification and visualisation 
df_pca <- df[, c("NuclRoundness")]

df_pca <- as.data.frame(df_pca)


## K-means clustering
set.seed(123) # always set specific seed in order to be able to replicate results

# Elbow approach
# Function to compute total within-cluster sum of square
wss <- function(k) {
  kmeans(df_pca, k, nstart = 10)$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15
# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares"
)
# determine perfect amount of clusters based on most extreme bend in plot 

# specify in centers how many clusters you have identified in the plot using the Elbow approach
km <- kmeans(df$NuclRoundness,centers=2)
df$clust_round <- as.factor(km$cluster)

# plot your data using the determined classification
ggplot(df, aes(x=NuclRoundness)) + 
  geom_histogram(aes(fill=clust_round),
                 binwidth=0.01)+
  stat_density(geom="line", color="red")

# in example - smaller value means nuclei are not perfectly round 
# two distributions visible 
# choose cluster with larger values 


##########
# you can also try to confirm the number of clusters that you determined using the Elbow approach by calculating the Silhouette score
# but this depends on the size of your data set (for me, my laptop sometimes ran out of memory trying this out)

# Silhouette approach
# Function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(df_pca, centers = k, nstart = 25, iter.max = 10, algorithm="Lloyd")
  ss <- silhouette(km.res$cluster, dist(df_pca))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15
# Extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes"
)

fviz_nbclust(df_pca, kmeans, method = "silhouette")
# Determine highest silhouette score
##########



# based on your plots, determine the number of your target cluster based on which you want to make the selection 
### double check before running !!!!!!!!!
# as you re-run the script the numbering of your target cluster might change 

df_cluster2 <- subset(df, df$clust_area == 2 & df$clust_hoechst == 1 & df$clust_round ==1)



## optional: 4.) MAP2 intensity 
# if you performed a MAP2 staining as well 

# (not part of example data set)

# you can try to identify clusters using Hierarchical clustering 
## Hierarchical clustering
# Identification and visualisation 
df_pca <- df_cluster2[, c("MAP2_MeanInten")]

df_pca <- as.data.frame(df_pca)


## K-means clustering
set.seed(123) # always set specific seed in order to be able to replicate results

# Elbow approach
# Function to compute total within-cluster sum of square
wss <- function(k) {
  kmeans(df_pca, k, nstart = 10)$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15
# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares"
)
# determine perfect amount of clusters based on most extreme bend in plot 


# specify in centers how many clusters you have identified in the plot using the Elbow approach
km <- kmeans(df_cluster2$MAP2_MeanInten,centers=2)
df_cluster2$clust <- as.factor(km$cluster)

#plot histogram with determined MAP2 clusters
ggplot(df_cluster2, aes(x=MAP2_MeanInten)) + 
  geom_histogram(aes(fill=clust))+
  stat_density(geom="line", color="red")


# BUT when I was plotting my MAP2 variable, I always found a bimodal distribution with a lot of values around 0 (MAP2 negative cells)
# I use the minimum here to differentiate these two peaks 
# because the separation according to clusters hierarchical clustering was not satisfactory 
find_local_minima_df <- function(df, breaks = 200, span = 0.3) {
  if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
  if (!require("zoo")) install.packages("zoo", dependencies = TRUE)
  library(dplyr)
  library(zoo)
  
  numeric_vars <- df %>% select(where(is.numeric))
  minima_list <- list()
  
  for (var_name in names(numeric_vars)) {
    data <- numeric_vars[[var_name]]
    hist_data <- hist(data, breaks = breaks, plot = FALSE)
    
    counts <- hist_data$counts
    mids <- hist_data$mids
    
    # Smooth the counts using a moving average (or LOWESS/LOESS if desired)
    smoothed_counts <- zoo::rollmean(counts, k = 3, fill = NA, align = "center")
    
    # Find local minima (where a point is less than its neighbors)
    minima_indices <- which(
      !is.na(smoothed_counts) &
        smoothed_counts < dplyr::lag(smoothed_counts) &
        smoothed_counts < dplyr::lead(smoothed_counts)
    )
    
    minima_midpoints <- mids[minima_indices]
    
    minima_list[[var_name]] <- minima_midpoints
  }
  
  return(minima_list)
}


find_local_minima_df(df_cluster2)

# determine the minimium and use it as a cut off to subset your data 
df_cluster2 <- subset(df_cluster2, df_cluster2$MAP2_MeanInten > 275)


#### Determine GFP pos vs GFP neg cells

# Identification and visualisation
df_pca <- df_cluster2[, c("GFP_MeanInten")]

df_pca <- as.data.frame(df_pca)


## K-means clustering
set.seed(123) # always set specific seed in order to be able to replicate results

# Elbow approach
# Function to compute total within-cluster sum of square
wss <- function(k) {
  kmeans(df_pca, k, nstart = 10)$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15
# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares"
)
# determine perfect amount of clusters based on most extreme bend in plot 

# specify in centers how many clusters you have identified in the plot using the Elbow approach
# in example - chose 3 clusters as boxplot showed these really extreme values which are nicely seperated with three clusters 
# but you could also argue to only use 2 clusters based on Elbow approach 
# as always there are multiple possibilities, but you need to be able to justify your choice 
km <- kmeans(df_cluster2$GFP_MeanInten,centers=3)
df_cluster2$clust_gfp <- as.factor(km$cluster)

# plot your data using the determined classification
ggplot(df_cluster2, aes(x=GFP_MeanInten)) + 
  geom_histogram(aes(fill=clust_gfp), bins = 100)+
  stat_density(geom="line", color="red")

boxplot(df_cluster2$GFP_MeanInten)


## DOUBLE CHECK NUMBERS OF CLUSTERS!!! see comment above
# identify GFP populations 
df_cluster2$GFPstatus <- ifelse(df_cluster2$clust_gfp == 3, "pos", ifelse(df_cluster2$clust_gfp == 2, "neg", "pospos"))
df_cluster2$GFPstatus <- as.factor(df_cluster2$GFPstatus)


df_cluster2 <- na.omit(df_cluster2)


# plot clusters with your staining intensity of protein of interest 
# in this example I was interested in Drebrin (my favourite protein on this planet) intensity 
ggplot(df_cluster2, aes(x=Dbn1_MeanInten)) + 
  geom_histogram(aes(fill=GFPstatus), bins = 100)+
  stat_density(geom="line", color="red")


# store your processed data using a file name you like 
write.table(df_cluster2,
            file = "SINEUP_Dbn1_N2ACells_scDATA_processeddata.csv",
            sep = ";", dec = ",", row.names = F
)



## The data is now ready for statistical analysis
## Multilevel analyses

# Intercept only model 
Dbn1_0 <- glm(Dbn1_MeanInten ~ 1, data = df_cluster2, family = Gamma(link="log")) 
summary(Dbn1_0)

# Random Intercept model 
Dbn1_1 <- glmer(Dbn1_MeanInten ~ 1 + (1 |ScreenName/WellPlate_Well), data = df_cluster2, family = Gamma(link="log")) 
summary(Dbn1_1)
Dbn1_2 <- glmer(Dbn1_MeanInten ~ 1 + (1 |ScreenName), data = df_cluster2, family = Gamma(link="log")) 
summary(Dbn1_2) 
Dbn1_3 <- glmer(Dbn1_MeanInten ~ 1 + (1 |WellPlate_Well), data = df_cluster2, family = Gamma(link="log")) 
summary(Dbn1_3)  

anova(Dbn1_3, Dbn1_0, test="Chisq")

# Factor model 
Dbn1_4 <- glmer(Dbn1_MeanInten_corr ~ 1 + GFPstatus * Condition + (1 |WellPlate_Well), data = df_cluster2, family = Gamma(link="log")) 
summary(Dbn1_5)  

anova(Dbn1_4, Dbn1_3, test="Chisq") 









