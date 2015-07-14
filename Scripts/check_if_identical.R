
identical_ans<-function(ans1,ans2) {

for (i in c(1,3,4,5))
{
ans1<-ans1[sort(rownames(ans1)),]
ans2<-ans2[sort(rownames(ans2)),]

print(identical(ans1[,i],ans2[,i]))

}

print("p-values")
print(cbind(round(ans1[,2],2),round(ans2[,2],2)))
print("qvalues")

print(cbind(round(ans1[,6],2),round(ans2[,6],2)))
}



ans1<-read.csv("LUAD_Neigh_45_3_MUT.matrix.csv-174935-2015-06-19_results_final.csv",row.names=1)

ans2<-read.csv("LUAD_Neigh_45_3_Mut.matrix.csv-180908-2015-06-17_results_final.csv",row.names=1)

rows<-intersect(rownames(ans1),rownames(ans2))

identical_ans(ans1[rows,],ans2[rows,])


pii1<-fread("LUAD_Neigh_45_3_MUT.matrix.csv-174935-2015-06-19_pii_values.csv",data.table=FALSE)
pii1<-pii1[,-1]
pii2<-fread("LUAD_Neigh_45_3_Mut.matrix.csv-180908-2015-06-17_pii_values.csv",data.table=FALSE)
pii2<-pii2[,-1]

#
print("Check if pii_files are identical:")
identical(pii1[,rows],pii2[,rows])







rows[1]






