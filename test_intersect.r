TestIntersect<-function(intv1,intv2){
	if(length(intv1)!=2&&length(intv2)!=2){
		return(NULL);
	}
	if(intv1[2]<intv2[1]||intv2[2]<intv1[1]){
		return(FALSE);
	}
	else{
		return(TRUE);
	}
}

FindIntersect<-function(intvs,intv){
	intvs<-as.matrix(intvs);
	if(ncol(intvs)!=2){
		return(NULL);
	}
	return(which(apply(intvs,1,TestIntersect,intv2=intv)));
}	

LocateGenes<-function(gene,intvs){
	gene.ps<-gene[,2:3];
	return(FindIntersect(intvs,gene.ps));
}

ReadGene<-function(gene_file){
	gene_data<-read.table(gene_file,sep='\t',header=T,stringsAsFactor=F,quote='')[,c(1,3:5)];
	return(gene_data);	
}

MapGenes<-function(chr_id,region.start,region.end,stat_df,gene_data){
	sub_gene_data<-subset(gene_data,chr==chr_id)[,c(1,3:4)];
	sub_gene_data.pos<-sub_gene_data[,2:3];
	sub_gene_data.nam<-sub_gene_data[[1]];
	region.covergenes.id<-FindIntersect(sub_gene_data.pos,c(region.start,region.end));
	if(length(region.covergenes.id)==0){
		return(list());
	}
	region.covergenes<-sub_gene_data[region.covergenes.id,];
	gene.coverage<-list();
	for(i in 1:length(region.covergenes.id)){
		tmp<-LocateGenes(region.covergenes[i,],stat_df[,1:2]);
		if(length(tmp)==0){
			next;
		}
		gene.coverage[[region.covergenes[i,1]]]<-tmp;	
	}
	return(gene.coverage);
}
