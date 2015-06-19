
identical_ans<-function(ans1,ans2) {

for (i in c(1,3,4,5))
{
ans1<-ans1[sort(rownames(ans1)),]
ans2<-ans2[sort(rownames(ans2)),]

print(identical(ans1[,i],ans2[,i]))

}

print("p-values")
print(round(ans1[,2],2))
print(round(ans2[,2],2))
print("qvalues")
print(ans1[,6])
print(ans2[,6])
}
