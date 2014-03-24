########################################
########################################
##	we will use package ape to plot all
##	of the coalescent reconstructed trees
##  library "ape" provide function of plot
##  tree and library "shape" provide function
##  to generate colors and basic plot elements
library(ape)
library(shape)
source('plot_chrome.r')
source('test_intersect.r')
source('get_mid.r')
###############################################################
##	initial set of some plot parameters
pos.title<-c(0.05,0.95,0.91,0.98);##  where the title region located
pos.ygene<-c(0.78,0.89);		  ##  yaxis of genes
pos.ytree<-c(0.38,0.77);		  ##  yaxis of Trees
pos.ysign<-c(0.17,0.35);		  ##  yaxis of signals trend	
pos.chrome<-c(0.02,0.98,0.08,0.14);##  where the chromosome region

dir.trees<-'downwards';## direction of trees should be vertical
###############################################################
select.col<-function(vec.select,col.bg='blue',col.ed='red'){
	vec.num<-length(vec.select);
	col.num<-ifelse(vec.num>20,vec.num,20);
	col<-shadepalette(col.num+1,col.bg,col.ed);
	vec.col<-seq(min(vec.select),max(vec.select),length.out=col.num+1);
	find.closest<-function(m,arrays)return(which.min(abs(m-arrays)));
	vec.order<-sapply(vec.select,find.closest,vec.col);
	col.ret<-col[vec.order];
	return(col.ret);
}

generate.col<-function(dt,key,col.bg='blue',col.ed='red'){
	dt$col<-select.col(dt[[key]],col.bg,col.ed);
	return(dt);
}

gen.treeind<-function(dt,key,M=20){
	ind.max<-which.max(dt[[key]]);
	n<-length(dt[[key]]);
	if(n<=M){
		return(1:n);
	}
	indx<-c();
	if(ind.max==1){
		indx<-c(1,get.mid(2,n,M-1));
	}
	else if(ind.max==n){
		indx<-c(get.mid(1,n-1,M-1),n);
	}
	else{
		p1<-1:(ind.max-1);
		p2<-(ind.max+1):n;
		p1.len<-length(p1);
		p2.len<-length(p2);
		p<-p1.len+p2.len;
		if(p1.len>=p2.len){
			n1<-ceiling(p1.len*(M-1)/p);
			n2<-floor(p2.len*(M-1)/p);
		}
		else{
			n1<-floor(p1.len*(M-1)/p);
			n2<-ceiling(p2.len*(M-1)/p);
		}
		indx<-c(get.mid(1,ind.max-1,n1),ind.max,get.mid(ind.max+1,n,n2));
	}
	return(indx);
}
			
gen.treepos<-function(TreeUnit,i,pos.axisbegin,pos.treebegin){
	if(i==1){
		return(c(pos.axisbegin,pos.treebegin+TreeUnit));
	}	
	else{
		return(c(pos.treebegin+(i-1)*TreeUnit,pos.treebegin+i*TreeUnit));
	}	
}

gen.recpos<-function(rec.beg,rec.end,pos.treebegin,rec.ini,scale){
	rec.l<-pos.treebegin+(rec.beg-rec.ini)*scale;
	rec.r<-pos.treebegin+(rec.end-rec.ini)*scale;
	return(c(rec.l,rec.r));
}

gen.sigpos<-function(pos.sig,pos.treebegin,rec.ini,scale){
	return(pos.treebegin+(pos.sig-rec.ini)*scale);
}

plot.empty<-function(x=0,y=0,title=NA,xlim=c(-1,1),ylim=c(-1,1)){
	plot(x,y,type='n',axes=F,xlab=NA,ylab=NA,main=title,xlim=xlim,ylim=ylim);
}

plot.forest<-function(filename,region){
	print(filename);
	#print(region);
########################################################
########################################################
	len.treeintv<-0.85;
	len.treeaxis<-0.005;
	pos.axisbegin<-0.045;
	pos.treebegin<-pos.axisbegin+len.treeaxis;
	pos.axisend<-pos.treebegin+len.treeintv;
	pos.treeend<-pos.axisend;
	base.num<-14;## all the above parameters is under the setting of base.num when it equals 18
###############################################################
	tree.height<-region[[1]][['height']];
	tree.thre<-quantile(tree.height,0.9);
	ylim.tree<-min(6000,tree.thre);
###############################################################
#### preprocessing the overall prameters
	forest<-region[[2]];
	size<-length(forest);
	scale.tree<-size/base.num;
	if(scale.tree<1.0){
		len.treeintv<-len.treeintv*scale.tree;
		len.treeaxis<-len.treeaxis*scale.tree;
	
		pos.treebegin<-(0.9-len.treeintv)/2;
		pos.axisbegin<-pos.treebegin-len.treeaxis;
	}
	#ylim.tree<-max(region[[1]][['time']])*2.7;
	if(length(region[[4]])==0){##case for no gene included
		pos.ytree<-c(0.43,0.88);                
		pos.ysign<-c(0.20,0.42);
		pos.chrome<-c(0.02,0.98,0.08,0.17);
	}
########################################################
########################################################
#### begin plot
	pdf(paste(filename,'.pdf',sep=''));
########################################################
########################################################
#### this part plot the title	
	par(plt=pos.title,cex=0.8);
	plot.empty();
	text(0,0,'Trend of Trees',col='darkgreen',cex=1.5,font=4);
########################################################
########################################################
#### this region plot the phylogentic trees
	tree.index<-gen.treeind(region[[1]],'h2');
	t.size<-length(tree.index);
	delta=len.treeintv/t.size;
	##############################################################
	abs.treeid<-c();
	abs.treecoordfx<-c();
	abs.treecoordfy<-c();
	##############################################################
	for(i in seq(tree.index)){
		id<-tree.index[i];
		plt.tree<-c(gen.treepos(delta,i,pos.axisbegin,pos.treebegin),pos.ytree);
		if(i==1){
			par(new=T,plt=plt.tree,cex=0.8);
			plot.phylo(forest[[id]],direction=dir.trees,show.tip.label=F,y.lim=ylim.tree,edge.color=region[[1]][['col']][id]);
			axisPhylo(2,mgp=c(0,0.3,0));
		}
		else{
			par(new=T,plt=plt.tree,pos.ytree);
			plot.phylo(forest[[id]],direction=dir.trees,show.tip.label=F,y.lim=ylim.tree,edge.color=region[[1]][['col']][id]);
		}
		#############################################################
		abs.treeid<-c(abs.treeid,id);
		abs.treecoordfx<-c(abs.treecoordfx,mean(plt.tree[1:2]));
		abs.treecoordfy<-c(abs.treecoordfy,pos.ytree[1]);
		#############################################################
	}
	bar.l=pos.treebegin+len.treeintv+0.15*delta;
	bar.r=bar.l+delta*0.35;
	bar.b=pos.ytree[1]+(pos.ytree[2]-pos.ytree[1])*0.2;
	bar.t=pos.ytree[1]+(pos.ytree[2]-pos.ytree[1])*0.6;
	par(new=T,plt=c(bar.l,bar.r,bar.b,bar.t));
	image.ymin=min(region[[1]][['coef']]);
	image.ymax=max(region[[1]][['coef']]);
	eps=0.00001;
	if(abs(image.ymax-image.ymin)<eps){
		image.ymin=0.75*image.ymin;
		image.ymax=1.25*image.ymax;
	}
	image.y=seq(image.ymin,image.ymax,length.out=ifelse(size>20,size+1,21));
	image.ylen=length(image.y);
	image(y=image.y,z=matrix(image.y,nrow=1),col=shadepalette(length(image.y),'red','orange'),xaxt='n',yaxt='n',ylab=NA);
	mtext('s',side=3,col='darkblue',font=4,cex=0.7,line=0.3);
	axis(side=4,tck=-0.07,labels=NA);
	axis.at=round(image.y[c(3,floor(image.ylen/2),image.ylen)],3);
	axis(side=4,at=axis.at,labels=as.character(axis.at),lwd=0,line=-0.98,cex.axis=0.7);

########################################################
#### this part plot the trend of signals	
	signals.pos<-rowMeans(region[[1]][,1:2]);
	signals.beg<-region[[1]][['begin']];
	signals.end<-region[[1]][['end']];
	signals.h2<-region[[1]][['h2']];
	signals.h1<-region[[1]][['h1']];
	signals.len<-length(signals.h2);
	h2.color='blue';
	h1.color='green';
	h.cex=0.8;
	signals.xl<-region[[1]][['begin']][1];
	signals.xr<-region[[1]][['end']][signals.len];
	signals.range<-signals.xr-signals.xl;
	signals.scale<-len.treeintv/signals.range;
	ylim.l<-min(region[[1]][['h1']])-3;
	ylim.r<-max(region[[1]][['h2']])+3;
	rect.max<-max(region[[1]][['time']])+5;
	rect.min<-0;
	rect.ceil<-ylim.r-1;
	rect.floor<-ylim.l+1.5;
	plt.par<-c(pos.axisbegin,pos.treeend,pos.ysign);
	par(new=T,plt=plt.par,cex=h.cex,xaxs='i');
	
	for(i in 1:signals.len){
		recty<-c();
		if(i==1){
			plot(signals.pos[i],signals.h2[i],type='p',col=h2.color,xlab=NA,ylab=NA,cex=2*h.cex,pch=18,axes=F,srt=45,xlim=c(signals.xl-2,signals.xr+2),ylim=c(ylim.l,ylim.r));
			
			points(signals.pos[i],signals.h1[i],col=h1.color,cex=1.7*h.cex,pch=19);
			axis(side=1,mgp=c(0,-0.1,-0.6),tck=0.005,lwd=0);
			axis(side=2,at=round(seq(ylim.l,ylim.r,(ylim.r-ylim.l)/5)),mgp=c(0,0.3,0));
			usr.l<-par('usr');
			text(usr.l[1],usr.l[4],'h1/h2',adj=c(0,1),cex=0.8,font=3);
			rect.height<-(rect.ceil-rect.floor)*(region[[1]][['time']][i]-rect.min)/(rect.max-rect.min)+rect.floor;
			recty<-c(recty,rect.height);
			rect(signals.beg[i],rect.floor,signals.end[i],rect.height,density=30,col=region[[1]][['col']][i]);
			if(signals.len==1){
				rect.scale=(rect.max-rect.min)/(rect.ceil-rect.floor);
				time.axis.at=round(seq(ylim.l,ylim.r,(ylim.r-ylim.l)/4));
				time.axis.label=round((seq(ylim.l,ylim.r,(ylim.r-ylim.l)/4)-ylim.l)*rect.scale+rect.min);
				axis(side=4,tck=-0.03,labels=NA);
				axis(side=4,line=-1.08,cex.axis=0.73,at=time.axis.at,labels=time.axis.label,lwd=0);	
				usr.r<-par('usr');
				text(usr.r[2],usr.r[4],'time',adj=c(1,1),cex=0.8,font=3);
				bar.l=0.915;
				bar.r=0.985;
				bar.b=pos.ysign[1]+(pos.ysign[2]-pos.ysign[1])*0.4;
				bar.t=pos.ysign[1]+(pos.ysign[2]-pos.ysign[1])*0.8;
				lines(signals.pos,signals.h2,col=h2.color,lwd=1.4);
				lines(signals.pos,signals.h1,col=h1.color,lwd=1.4);
				######################################################
				abs.sigbeg<-grconvertX(signals.beg,'user','ndc');
				abs.sigend<-grconvertX(signals.end,'user','ndc');
				abs.treecoordtx<-grconvertX(signals.pos,'user','ndc');
				abs.sigpy1<-grconvertY(signals.h2,'user','ndc');
				abs.sigpy2<-grconvertY(recty,'user','ndc');
				abs.treecoordty<-(abs(abs.sigpy1+abs.sigpy2)+abs(abs.sigpy1-abs.sigpy2))/2;
				######################################################
				par(new=T,plt=c(bar.l,bar.r,bar.b,bar.t),xaxs='i');
				plot(-0.2,0.8,xlab=NA,ylab=NA,axes=F,cex=1.5*h.cex,pch=18,col=h2.color,xlim=c(-0.52,0.92),ylim=c(-1,1));
				text(0.3,0.8,'h2',col=h2.color,cex=h.cex);
				points(-0.2,0.4,cex=1.2*h.cex,pch=19,col=h1.color);
				segments(-0.23,0.8,-0.18,0.8,col=h2.color,lwd=2.3);
				segments(-0.23,0.4,-0.18,0.4,col=h1.color,lwd=2.3);
				text(0.3,0.4,'h1',col=h1.color,cex=h.cex);
				rect(-0.35,-0.8,-0.10,0.0,density=30,col='red');
				text(0.35,-0.45,'time',col='red',cex=1.1*h.cex);	
			}
		}
		else{
			points(signals.pos[i],signals.h2[i],col=h2.color,cex=2*h.cex,pch=18);
			points(signals.pos[i],signals.h1[i],col=h1.color,cex=1.7*h.cex,pch=19);
			rect.height<-(rect.ceil-rect.floor)*(region[[1]][['time']][i]-rect.min)/(rect.max-rect.min)+rect.floor;
			recty<-c(recty,rect.height);
			if(i!=signals.len){
				rect(signals.beg[i],rect.floor,signals.end[i],rect.height,density=30,col=region[[1]][['col']][i]);
			}
			else{
				rect(signals.beg[i],rect.floor,signals.end[i],rect.height,density=30,col=region[[1]][['col']][i]);
				rect.scale=(rect.max-rect.min)/(rect.ceil-rect.floor);
				time.axis.at=round(seq(ylim.l,ylim.r,(ylim.r-ylim.l)/4));
				time.axis.label=round((seq(ylim.l,ylim.r,(ylim.r-ylim.l)/4)-ylim.l)*rect.scale+rect.min);
				axis(side=4,tck=-0.03,labels=NA);
				axis(side=4,line=-1.08,cex.axis=0.73,at=time.axis.at,labels=time.axis.label,lwd=0);	
				usr.r<-par('usr');
				text(usr.r[2],usr.r[4],'time',adj=c(1,1),cex=0.8,font=3);
				bar.l=0.915;
				bar.r=0.985;
				bar.b=pos.ysign[1]+(pos.ysign[2]-pos.ysign[1])*0.4;
				bar.t=pos.ysign[1]+(pos.ysign[2]-pos.ysign[1])*0.8;
				lines(signals.pos,signals.h2,col=h2.color,lwd=1.4);
				lines(signals.pos,signals.h1,col=h1.color,lwd=1.4);
				######################################################
				abs.sigbeg<-grconvertX(signals.beg,'user','ndc');
				abs.sigend<-grconvertX(signals.end,'user','ndc');
				abs.treecoordtx<-grconvertX(signals.pos,'user','ndc');
				abs.sigpy1<-grconvertY(signals.h2,'user','ndc');
				abs.sigpy2<-grconvertY(recty,'user','ndc');
				abs.treecoordty<-(abs(abs.sigpy1+abs.sigpy2)+abs(abs.sigpy1-abs.sigpy2))/2;
				######################################################
				par(new=T,plt=c(bar.l,bar.r,bar.b,bar.t),xaxs='i');
				plot(-0.2,0.8,xlab=NA,ylab=NA,axes=F,cex=1.5*h.cex,pch=18,col=h2.color,xlim=c(-0.52,0.92),ylim=c(-1,1));
				text(0.3,0.8,'h2',col=h2.color,cex=h.cex);
				points(-0.2,0.4,cex=1.2*h.cex,pch=19,col=h1.color);
				segments(-0.23,0.8,-0.18,0.8,col=h2.color,lwd=2.3);
				segments(-0.23,0.4,-0.18,0.4,col=h1.color,lwd=2.3);
				text(0.3,0.4,'h1',col=h1.color,cex=h.cex);
				rect(-0.35,-0.8,-0.10,0.0,density=30,col='red');
				text(0.35,-0.45,'time',col='red',cex=1.1*h.cex);	
			}
		}	
	}	
########################################################
########################################################
#### this part plot the genes location if there are
	abs.treecoordtx<-abs.treecoordtx[abs.treeid];
	abs.treecoordty<-abs.treecoordty[abs.treeid];
	if(length(region[[4]])!=0){
		genes<-region[[4]];
		gene.names<-names(genes);
		gene.len<-length(gene.names);
		gene.bot<-pos.ygene[1];
		gene.top<-pos.ygene[2];
		gene.height<-seq(gene.bot,gene.top,length.out=gene.len+2);
		gene.height<-gene.height[c(-1,-length(gene.height))];
		gene.ydelta<-0.33*(gene.top-gene.bot)/(gene.len+1);
		gene.ysub<-cbind(gene.height-gene.ydelta,gene.height+gene.ydelta);
		for(i in 1:nrow(gene.ysub)){
			genes.id<-genes[[gene.names[i]]];
			genes.num<-length(genes.id);
			genes.col<-region[[1]][['col']][genes.id];
			genes.lpos<-abs.sigbeg[genes.id[1]];
			genes.rpos<-abs.sigend[genes.id[genes.num]];
			par(new=T,plt=c(genes.lpos,genes.rpos,gene.ysub[i,]),xaxs='i');
			plot(0,0,xlab=NA,ylab=NA,xlim=c(-1,1),ylim=c(-0.5,0.5),axes=F,type='n');
			textflag(mid=c(0,0),radx=0.99,rady=0.4,rx=0.01,lab=gene.names[i],cex=0.7,font=2,col=genes.col);	
		}		
	}
#############################################################
#############################################################
#### add arrows and connect component in different layers
	par(new=T,plt=c(0,1,0,1),xaxs='i',yaxs='i');
	plot.empty(xlim=c(0,1),ylim=c(0,1));
	arrow.len<-length(abs.treeid);
	#for(i in 1:arrow.len){
		#arrows(abs.treecoordfx[i],abs.treecoordfy[i],abs.treecoordtx[i],abs.treecoordty[i],col=region[[1]][['col']][i],lwd=0.8);
		Arrows(abs.treecoordfx,abs.treecoordfy,abs.treecoordtx,abs.treecoordty,lwd=1.2,col=region[[1]][['col']][abs.treeid],arr.adj=1,arr.length=0.2,arr.width=0.05,arr.type='simple');
	#}
#############################################################
#### this region plot the chromosome	
	par(new=T,plt=pos.chrome,xaxs='i');
	range.bg=(pos.treebegin-pos.chrome[1])/(pos.chrome[2]-pos.chrome[1]);
	range.ed=(pos.treebegin+len.treeintv-pos.chrome[1])/(pos.chrome[2]-pos.chrome[1]);
	plot.chrome(range.bg,range.ed,region[[3]][1],loc=region[[3]][2:3],loc.col='red');
########################################################
########################################################	
	dev.off();
}

get.par<-function(dt,ResFormats,TreeFormats,IntFormats,gene.data){
	region.chr<-dt[['chr']][1];
	region.start<-dt[['start']][1];
	region.end<-dt[['end']][nrow(dt)];
	region.posbg<-dt[['begin']];
	region.posed<-dt[['end']];
	region.estim<-dt[['est.time']];
	region.escoe<-apply(cbind(dt[['exp.lm.e']],dt[['det.lm.d']]),1,min);
	region.bstrs<-paste('\\[',region.posbg,sep='');
	region.index<-sapply(region.bstrs,function(m,tf) return(grep(m,tf[[1]])),tf=TreeFormats);
	region.treef<-TreeFormats[region.index,2];
	region.resut<-ResFormats[region.index,c(2,3,7)];
	int.len<-ncol(IntFormats);
	region.height<-apply(IntFormats[region.index,2:int.len],1,sum);
	options(stringsAsFactors=F);
	region.dt<-data.frame(region.posbg,region.posed,region.estim,region.escoe,region.treef,region.resut[,1:3],region.height);
#########################################################
#	region.dt<-region.dt[!duplicated(region.dt[3]),];
	eliminate_dup<-function(region.dt){
		region.dt.stack=region.dt[1,];
		i<-1;
		dt.count<-1;
		while(i<nrow(region.dt)){
			j=i+1;
			while(region.dt[i,5]==region.dt[j,5]){
				if(j==nrow(region.dt)){
					region.dt.stack[dt.count,2]=region.dt[j,2];
					return(region.dt.stack);
				}
				else{
					j=j+1;
				}
			}
			if(j!=i+1){
				region.dt.stack[dt.count,2]=region.dt[j-1,2];
			}
			region.dt.stack=rbind(region.dt.stack,region.dt[j,]);
			dt.count=dt.count+1;
			i=j
		}
		return(region.dt.stack);
	}		

	region.dt<-eliminate_dup(region.dt);
#########################################################
	region.cdt<-as.character(region.dt[[5]]);
	region.trees<-lapply(region.cdt,function(m) return(read.tree(text=m)));	
	region.dt[[5]]<-NULL;
	colnames(region.dt)<-c('begin','end','time','coef','a','h1','h2','height');
	region.dt<-generate.col(region.dt,'coef',col.bg='red',col.ed='orange');
	region.genes<-MapGenes(region.chr,region.start,region.end,region.dt,gene.data);	
	region.list<-list(info=region.dt,tree=region.trees,bs=c(region.chr,region.start,region.end),genes=region.genes);
	return(region.list);
}

makeplot<-function(CandGene,Res.files,Tree.files,Int.files,prefix,gene.data,output.dir){
	mdt<-CandGene;
	colnames(mdt)[3:4]<-c('begin','end');
	match1.file<-function(m,f) return(f[grep(paste('chr',m,'_',sep=''),f)]);
	match2.file<-function(m,f)return(f[grep(paste('chr',m,'\\.',sep=''),f)]);
	res.file<-match1.file(mdt[[1]][1],Res.files);
	tree.file<-match2.file(mdt[[1]][1],Tree.files);
	int.file<-match1.file(mdt[[1]][1],Int.files);
	resf<-read.table(res.file);
	treef<-read.table(tree.file,stringsAsFactors=F);
	intf<-read.table(int.file);
	tp<-get.par(mdt,resf,treef,intf,gene.data);
	filename<-paste(output.dir,'chr',mdt[[1]][1],'_',mdt[[2]][1],sep='');
	plot.forest(filename,tp);
}

plot.pop<-function(input.dir,prefix,gene_file,output.dir){
	input.dir<-gsub('\\/$','',input.dir);
	output.dir<-gsub('(\\/)?$','/',output.dir);
	dirs<-list.dirs(input.dir,full.names=T,recursive=F);
	dirs<-dirs[grep(prefix,dirs)];
	file<-list.files(input.dir,full.names=T,pattern=paste('cand_',prefix,sep=''));
	gene<-ReadGene(gene_file);
	Cand<-read.table(file,header=T);
	res.dirs<-paste(dirs,'/','res',sep='');
	tree.dirs<-paste(dirs,'/','Tree',sep='');
	int.dirs<-paste(dirs,'/','tree',sep='');
	res.files<-list.files(res.dirs,full.names=T);
	tree.files<-list.files(tree.dirs,full.names=T);
	int.files<-list.files(int.dirs,full.names=T);
	by(Cand,list(Cand$chr,Cand$start),makeplot,Res.files=res.files,Tree.files=tree.files,Int.files=int.files,prefix=prefix,output.dir=output.dir,gene.data=gene);
}

args<-commandArgs(T);
indir<-args[1];
prefix<-args[2];
gene_file<-args[3];
outdir<-args[4];
plot.pop(indir,prefix,gene_file,outdir);
q(save='no');
