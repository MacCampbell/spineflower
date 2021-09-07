### Script adapted from Matteo Fumagalli ###

args=commandArgs(T)
fin=args[1]
nPop1=as.numeric(args[2])*2
nPop2=as.numeric(args[3])*2
pop1=args[4];
pop2=args[5];


rm(args)

ANGSD.2D.SFS <- scan(paste(fin, sep=""), quiet=T)
ANGSD.2D.SFS <- t(matrix(ANGSD.2D.SFS, nrow=nPop2+1, ncol=nPop1+1))

n1=nrow(ANGSD.2D.SFS)-1;
n2=ncol(ANGSD.2D.SFS)-1;

colv<-rainbow(n=40,start=0, end =0.3)

fout=paste(fin,".2D.pdf", sep="", collapse="")

pdf(file=fout)
image(x=seq(0,n1), y=seq(0,n2), col=colv, z=-log10(ANGSD.2D.SFS), xlab=pop1, ylab=pop2, main="Joint SFS")
dev.off()

cat("Output file:", fout, "\n")

