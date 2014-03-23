#***********************************************************************
#CHR INFO DATA**********************************************************
Rchr.len<-scan("chr_len.txt",sep=" ");                                #**
scalef<-1.0/Rchr.len[1];											 #**
chr.len<-Rchr.len*scalef;											 #**
chr.cen<-read.table("genome_info.txt");								 #**
chr.cen<-chr.cen[,3:4]*scalef;                                   #**
chr.cen.mean<-(chr.cen[[1]]+chr.cen[[2]])/2.0;                       #**
chr.cen[[1]]<-chr.cen.mean-0.0001;chr.cen[[2]]<-chr.cen.mean+0.0001;	 #**
#***********************************************************************
#***********************************************************************
#Using roundrect function in package shape
library(shape)
#
#***********************************************************************
#***********************************************************************
get.scale.factor<-function(){
	return(scalef);
}

plot.chrome<-function(range.bg,range.ed,chr.id,mid.pos=0,col='lightblue',xlim=c(-1,1),ylim=c(-0.025,0.025),rady=0.01,mid.y=0,lim.scale=T,loc=NULL,loc.col='red',arrow.len=0.2,arrow.lwd=2.5,plot.arrow=T,lcol='pink'){
	chrome.len<-chr.len[chr.id];
	if(lim.scale) xlim=c(-chrome.len/2-0.05,chrome.len/2+0.05);
	chrome.cen.l<-chr.cen[[1]][chr.id];
	chrome.cen.r<-chr.cen[[2]][chr.id];
	chrome.b<-mid.pos-chrome.len/2;
	chrome.e<-mid.pos+chrome.len/2;
	chrome.cen.l<-chrome.b+chrome.cen.l;
	chrome.cen.r<-chrome.b+chrome.cen.r;
	mid.lx<-mean(c(chrome.b,chrome.cen.l));
	mid.rx<-mean(c(chrome.cen.r,chrome.e));
	radx.l<-(chrome.cen.l-chrome.b)/2-rady;
	radx.r<-(chrome.e-chrome.cen.r)/2-rady;
	plot(mid.pos,mid.y,xlab=NA,ylab=NA,axes=F,type='n',xlim=xlim,ylim=ylim);
	roundrect(c(mid.lx,mid.y),radx=radx.l,rady=rady,col=col,lcol=lcol);
	roundrect(c(mid.rx,mid.y),radx=radx.r,rady=rady,col=col,lcol=lcol);
	text(chrome.b,-rady-0.002,as.character(chr.id),font=2);
	if(!is.null(loc)){
		l.bg<-chrome.b+loc[1]*scalef;
		l.ed<-chrome.b+loc[2]*scalef;
		rect(l.bg,-rady,l.ed,rady,col=loc.col);
		if(plot.arrow){
			#arrow.xl<-seq(xlim[1],xlim[2],length.out=11)[4];
			#arrow.xr<-seq(xlim[1],xlim[2],length.out=11)[7];
			#arrow.xl<-xlim[1];
			#arrow.xr<-xlim[2];
			arrow.xl<-(xlim[2]-xlim[1])*range.bg+xlim[1];
			arrow.xr<-(xlim[2]-xlim[1])*range.ed+xlim[1];
			Arrows(arrow.xl,ylim[2],l.bg,rady+0.002,arr.length=arrow.len,arr.type='triangle',arr.col=loc.col,lwd=arrow.lwd);
			Arrows(arrow.xr,ylim[2],l.ed,rady+0.002,arr.length=arrow.len,arr.type='triangle',arr.col=loc.col,lwd=arrow.lwd);
	
		}
	}	
}
