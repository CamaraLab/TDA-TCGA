x<-read.csv("scan_summary.csv")


#Number of samples per graph plot:
ggplot(x, aes(x=factor(resolution), y=factor(gain), label=q_0.15,fill=q_0.15)) + 
  scale_fill_gradient2(low = 'maroon',high = 'blue') +
  geom_tile() + theme_minimal() + geom_text() + ggtitle("Genes q<=0.15")

y<-read.csv("scan_summary.csv")
ggplot(y, aes(x=factor(resolution), y=factor(gain), fill=mutload)) + 
  scale_fill_continuous(low = 'red',high = 'blue',guide_legend(title = "p_value")) +
  geom_tile() + theme_minimal() + geom_text(label=round(y$mutload,2),size=7,vjust=-0.1) + ggtitle("Mutational load - Colon") + 
  xlab("Resolution") + ylab ("Gain") + geom_text(label=paste0("(",y$first_connected_samples,")"),vjust=1,size=6) +
  coord_equal()
  

qplot(data=x,x=resolution,y=gain,geom="tile",fill=q_0.15,label=q_0.15)
x$gain

