x<-read.csv("scan_summary.csv")


#Number of samples per graph plot:
ggplot(x, aes(x=factor(resolution), y=factor(gain), label=q_0.15,fill=q_0.15)) + 
  scale_fill_gradient2(low = 'maroon',high = 'blue') +
  geom_tile() + theme_minimal() + geom_text() + ggtitle("Genes q<=0.15")

y<-read.csv("scan_summary.csv")
ggplot(y, aes(x=factor(resolution), y=factor(gain), label=round(mutload,2),color=mutload)) + 
  scale_color_gradient(low = 'red',high = 'blue',guide_legend(title = "p_value")) +
  geom_point(size=6) + theme_bw() + geom_text(vjust=1.6) + ggtitle("Mutational load - Colon") + 
  xlab("Resolution") + ylab ("Gain")
  

qplot(data=x,x=resolution,y=gain,geom="tile",fill=q_0.15,label=q_0.15)
x$gain

