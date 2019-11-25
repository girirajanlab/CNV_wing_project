#This R script performs k-means clustering of Drosophila adult wing phenotypic data for homologs of CNV genes and neurodevelopmental genes.
#Note: Code adapted from https://rpubs.com/williamsurles/310847

data<-read.table("wing_clustering_data_f.txt",sep="\t",header=T,row.names=1)
values<-subset(data,select= EV_0:BR_4)

#Determine best number of clusters, out of 1-15
wss <- 0
for (i in 1:15) {
     kmeans_cluster <- kmeans(values, centers = i, nstart = 20)
     wss[i] <- kmeans_cluster$tot.withinss
 }

#Plot results for all k-values and save as PDF. Pick a value that is at the "elbow" of the plot for # of clusters
plot(1:15, wss, type = "b", 
      xlab = "Number of Clusters", 
      ylab = "Within groups sum of squares")
#For this data, selected 5 as value based on scoring (no phenotype, mild, moderate, severe, lethal)

#Perform clustering of data using selected k-value (nstart=number of iterations for assigning values to cluster)
kmeans_cluster <- kmeans(values, centers = 5, nstart = 20)

#Get summary statistics to save in TXT file
kmeans_cluster

#Write cluster assignment to file
write.table(kmeans_cluster$cluster,file="kmeans_5_clusters.txt",quote=FALSE,row.names=TRUE,col.names=FALSE,sep="\t")

