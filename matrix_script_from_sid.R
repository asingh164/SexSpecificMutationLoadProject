                      ### Making the M matrix, inefficient solution... ###
reduced <- read.table("reduced.2L.txt", stringsAsFactors = F, header=T)

allele.matrix <- reduced[,c(5:354)] # This matrix should only contain allele data (should only have bases)

similarity.matrix <- matrix(nrow=ncol(allele.matrix),ncol=ncol(allele.matrix))
for(col in 1:ncol(allele.matrix)){
  matches <- allele.matrix[,col]==allele.matrix
  match.counts <- colMeans(matches, na.rm = T)
  similarity.matrix[,col] <- match.counts
}
colnames(similarity.matrix) <- colnames(allele.matrix)
rownames(similarity.matrix) <- colnames(allele.matrix)

write.table(similarity.matrix, file = 'M.matrix.reducedpk.2L.txt', quote=F, row.names=T, col.names=T, sep='\t')
                                   

				
				 ### Making the A matrix ###
reduced <- read.table("reducedpk.all.txt", stringsAsFactors = F, header=T)

ref <- reduced$V3
allele.matrix <- reduced[,5:354] # make sure this matrix contains only columns with allele data (should only have bases)
similarity.matrix <- matrix(nrow=ncol(allele.matrix),ncol=ncol(allele.matrix))
colnames(similarity.matrix) <- colnames(allele.matrix)
rownames(similarity.matrix) <- colnames(allele.matrix)


# Make pk variable
pk <- numeric(nrow(reduced))
matches <- ref==allele.matrix
pk <- rowMeans(matches, na.rm=T)

# Make Xij/Xji variables then cbind
rm(allele.matrix)
x.array <- matches*1

# Use -1, 1 instead of 0, 1.
  #x.array[x.array==0] <- -1

rm(matches)
x.array <- cbind(x.array, pk)



options(scipen=999)
# For loop to calculate matrix values
for (i in 1:ncol(similarity.matrix)){
  for (j in i:ncol(similarity.matrix)){
    similarity.matrix[i,j] <- similarity.matrix[j,i] <- 
    #A0 from paper: 
    mean((x.array[,i]-2*pk)*(x.array[,j]-2*pk)/(2*pk*(1-pk)), na.rm=T)
    #A0
    #mean((x.array[,i]-pk)*(x.array[,j]-pk)/(pk*(1-pk)), na.rm=T)  
    A
    sum((x.array[,i]-pk)*(x.array[,j]-pk), na.rm=T)/(sum(pk*(1-pk), na.rm=T))
  }
}
write.table(similarity.matrix, file = 'A0.paper_eqn.all.txt', quote=F, row.names=T, col.names=T, sep='\t')
