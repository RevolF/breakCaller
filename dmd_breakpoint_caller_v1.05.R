	
# newly added plot starts here, use histogram and set breaks independantly
# set span up and down to 15% areas around suspicious breakpoint

mainExc=function(filename,mark){
	
	library(Cairo)
	
	# avoid columns more than colnames by using header=F and fill=T
	rawDf=read.table(filename,stringsAsFactors=F,header=F,fill=T)
	
	splt=strsplit(filename,'_')[[1]]
	rawname=paste(splt[1],splt[2],sep='_')
	depthfile=paste(rawname,'_rawlib.basecaller.sorted.bam.dmddepth',sep='')
	
	depthdf=read.table(depthfile,stringsAsFactors=F,header=F)
	maxdepth=max(depthdf$V3,na.rm=T)
	
	if(!dir.exists(paste('./',rawname,sep=''))){
		dir.create(paste('./',rawname,sep=''))
		}
	if(!dir.exists(paste('./',rawname,'/',mark,sep=''))){
		dir.create(paste('./',rawname,'/',mark,sep=''))
		}
	
	for(i in 2:200){
		# print(rawDf$V5[i])
		tryCatch(
		{
		span=strsplit(rawDf$V1[i],'__')[[1]]
		
		if(rawDf$V2[i]>10 & span[1]!=span[2] & span[1]!='NA' & span[2]!='NA'){
		
			slide=round(0.15*as.numeric(rawDf$V5[i]))
			
			rgt=as.numeric(rawDf$V3[i])
			lft=as.numeric(rawDf$V4[i])
			
			if(as.numeric(rawDf$V3[i])<as.numeric(rawDf$V4[i])){
				lft=as.numeric(rawDf$V3[i])
				rgt=as.numeric(rawDf$V4[i])
				}
			
			tmpdpt=subset(depthdf,V2>(lft-slide) & V2<(rgt+slide))
			
			plotname=paste(rawDf$V2[i],rawDf$V1[i],rawDf$V3[i],rawDf$V4[i],rawDf$V5[i],mark,sep='_')
			
			CairoPNG(paste('./',rawname,'/',mark,'/',plotname,'.png',sep=''),width=1800,height=600)
			
			refbreaks=150
			# if(as.numeric(rawDf$V5[i])<150){
				# refbreaks=as.numeric(rawDf$V5[i])
				# }
			
			gap=round((abs(lft-rgt)+slide*2)/refbreaks)
			
			mat=matrix(ncol=4,nrow=refbreaks)
			
			for(k in 1:nrow(mat)){
				mat[k,1]=lft-slide+(k-1)*gap
				mat[k,2]=lft-slide+k*gap
				mat[k,3]=sum(as.numeric(subset(tmpdpt,V2>=mat[k,1] & V2<mat[k,2])[,3]),na.rm=T)
				}
			
			mat[,4]=ifelse(mat[,1]>lft & mat[,2]<rgt,5,8)
			
			barplot(beside=T,mat[,3],names.arg=mat[,1],col=mat[,4],main=plotname,xlab='chrX pos',ylab='depth',cex.lab=2.5,cex.main=2.5)
			
			dev.off()
			}
		},error=function(e){cat()}
		)
		}
	}

# Args=commandArgs()

# Args=commandArgs()
# filenamer=Args[4]
# markr=Args[5]
# mainExc(filename=filenamer,mark=markr)


tryCatch(
	{
	Args=commandArgs()
	filenamer=Args[4]
	markr=Args[5]
	mainExc(filename=filenamer,mark=markr)
	},error=function(e){cat()}
	)

	
	
	
	
	
	
	
	