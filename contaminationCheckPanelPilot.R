args = commandArgs(trailingOnly=TRUE)

library(matrixStats)
fileName=args[1]
print(fileName)
df = read.table(fileName,header=T,sep="\t")
first = T;
for(i in unique(df$chr)){
    if(first){
        res=diff(df[df$chr==i,]$tumorModifiedBAF);
        first=F;
    }
    else{
       res = c(res,diff(df[df$chr==i,]$tumorModifiedBAF))
    }
}

print( sd(res)/sqrt(2))
# proposed cut-off value = 0.05
