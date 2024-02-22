ORF.finder=function(sequences,target=NULL,expected.length=NULL,plus=0,minus=0,ANTISENS=F){
  T0=Sys.time()
  
  require(seqinr)
  require(dplyr)
  
  if(!is.null(expected.length)){Range=c((expected.length-minus),(expected.length+plus))}else{Range=NULL}
  if(ANTISENS){frames=c(0,1,2,0,1,2)}else{frames=c(0,1,2)}
  
  for (i in 1:length(sequences)){
    
    translation.table=NULL
    for (j in 1:length(frames)){
      if(j<4){S='F'}else{S='R'}
      AA=seqinr::translate(sequences[[i]], frame = frames[j], sens = S)
      
      
      if(AA[length(AA)]=='*'){AA=AA[1:(length(AA)-1)]}
      
      tmp=data.frame(
        Frame=frames[j],
        length=length(AA),
        Sens=S,
        Stop=ifelse('*' %in% AA,"STOP","GO")
      )
      if(j==1){translation.table=tmp}else{translation.table=rbind(translation.table,tmp)}
    }

    Best.AA=translation.table[translation.table$Stop=='GO',]

    if(nrow(Best.AA)==0){write.fasta(sequences = sequences[i], names = names(sequences)[i], file.out =  paste(target,"Non ORF DNA.fasta",sep=' '), open="a",nbchar = length(sequences[[i]]))}

    if(nrow(Best.AA)>=1){
      
      for(n in 1:nrow(Best.AA)){
        AA=seqinr::translate(sequences[[i]], frame = Best.AA$Frame[n], sens = Best.AA$Sens[n])
        if(!is.null(expected.length)){
          if(between(length(AA),min(Range),max(Range))){
            write.fasta(sequences = AA, names = paste(names(sequences)[i],paste('ORF',n,sep='_')), file.out =  paste(target,"ORF AA expected length.fasta",sep=' '), open="a",nbchar = length(AA))
            write.fasta(sequences = toupper(sequences[[i]]), names = names(sequences)[i], file.out =  paste(target,"ORF DNA expected length.fasta",sep=' '), open="a",nbchar = length(sequences[[i]]))
          }else{
            write.fasta(sequences = AA, names = paste(names(sequences)[i],paste('ORF',n,sep='_')), file.out =  paste(target,"ORF AA non expected length.fasta",sep=' '), open="a",nbchar = length(AA))
            write.fasta(sequences = toupper(sequences[[i]]), names = names(sequences)[i], file.out =  paste(target,"ORF DNA non expected length.fasta",sep=' '), open="a",nbchar = length(sequences[[i]]))
          }
        }else{
          write.fasta(sequences = AA, names = paste(names(sequences)[i],paste('ORF',n,sep='_')), file.out =  paste(target,"ORF AA.fasta",sep=' '), open="a",nbchar = length(AA))
          write.fasta(sequences = toupper(sequences[[i]]), names = names(sequences)[i], file.out =  paste(target,"ORF DNA.fasta",sep=' '), open="a",nbchar = length(sequences[[i]]))
        }
      }
    }
    
    Checkpoints=as.data.frame(quantile(1:length(sequences)))
    
    if(i==round(Checkpoints[,1][2],0)){print(paste('25% sequences processed; Time ellapsed:',round(as.numeric(Sys.time()-T0),2),'secs'))}
    if(i==round(Checkpoints[,1][3],0)){print(paste('50% sequences processed; Time ellapsed:',round(as.numeric(Sys.time()-T0),2),'secs'))}
    if(i==round(Checkpoints[,1][4],0)){print(paste('75% sequences processed; Time ellapsed:',round(as.numeric(Sys.time()-T0),2),'secs'))}
    if(i==round(Checkpoints[,1][5],0)){print(paste('100% sequences processed; Time ellapsed:',round(as.numeric(Sys.time()-T0),2),'secs'))}
    
  }
  
  if(!is.null(expected.length)){
    prot.file=seqinr::read.fasta(paste(target,"ORF AA expected length.fasta",sep=' '),seqtype = 'AA')
  }else{prot.file=seqinr::read.fasta(paste(target,"ORF AA.fasta",sep=' '))}
  
  ORFs=c()
  for(i in 1:length(prot.file)){ORFs[i]=unlist(strsplit(attr(prot.file[[i]],"Annot"),split=' '))[2]}
  
  if('ORF_2'%in%ORFs|'ORF_3'%in%ORFs|'ORF_4'%in%ORFs|'ORF_5'%in%ORFs|'ORF_6'%in%ORFs){
    warning(paste('Multiple ORFs detected for the following sequences: ',paste(names(prot.file)[which(ORFs!='ORF_1')],collapse='; '),sep="\n" ) ) }
  
}
check.BLASTp=function(AA.BLASTed,searchFor,To_exclude=NULL,DNA.BLASTed=NULL,nHits=1,DIR='.'){
  
  require(seqinr)
  
  if(is.character(DIR)){files=list.files(DIR);files=files[grep('-Alignment.txt',files)]}
  if(length(files)>0){

    if(length(files)>1){
      collated.files=NULL
      for (f in 1:length(files)){
        file=read.delim(paste(DIR,files[f],sep='/'),h=F)
        file=file[-grep('ASV_',file$V1)[1],'V1',drop=F]
        if(f==1){collated.files=file}else{collated.files=rbind(collated.files,file)}
      }
    }else{
      collated.files=read.delim(paste(DIR,files[1],sep='/'),h=F)
    }
    
    BLASTs=grep('ASV_',collated.files$V1)
    
    for (i in 1:length(BLASTs)){
      
      if(collated.files$V1[BLASTs[i]+1]!="No significant similarity found."){
        
        MATCH=F
        index=1
        while(MATCH==F){
          if(unlist(strsplit(collated.files$V1[BLASTs[i]+index],split=' '))[1]=='Description'){MATCH=T}
          if(MATCH==F){index=index+1}
             }

        for (n in 1:nHits){
          
          Hit=collated.files$V1[BLASTs[i]+(index+n)]
          match=F
          for(s in 1:length(searchFor)){
            if(identical(grep(searchFor[s],Hit),integer(0))==F){
              match=T
              break
            }
          }

          if(!is.null(To_exclude)){ 
            Hit2=gsub("\\[","",(gsub("\\]","",unlist(strsplit(collated.files$V1[BLASTs[i]+(3+n)],split=' ')))))
            if(identical(grep(TRUE,To_exclude%in%Hit2),integer(0))==F){match=F}
          }

          if(match){
            to_recover.AA=paste(unlist(strsplit(collated.files$V1[BLASTs[i]],split=' '))[3:4],collapse=' ')
            to_recover.DNA=unlist(strsplit(collated.files$V1[BLASTs[i]],split=' '))[3]
            seqinr::write.fasta(sequences =toupper(AA.BLASTed[[grep(to_recover.AA,names(AA.BLASTed))]]), 
                                names = names(AA.BLASTed)[grep(to_recover.AA,names(AA.BLASTed))], 
                                file.out =  paste(Target,"recovered ORF AA.fasta",sep=' '), open="a",nbchar = length(AA.BLASTed[[grep(to_recover.AA,names(AA.BLASTed))]]))
            
            if(!is.null(DNA.BLASTed)){
              seqinr::write.fasta(sequences =toupper(DNA.BLASTed[[  which(names(DNA.BLASTed)==to_recover.DNA)  ]]), 
                                  names = names(DNA.BLASTed)[  which(names(DNA.BLASTed)==to_recover.DNA) ], 
                                  file.out =  paste(Target,"recovered ORF DNA.fasta",sep=' '), open="a",
                                  nbchar = length(DNA.BLASTed[[  which(names(DNA.BLASTed)==to_recover.DNA)  ]]))}
            
            break
          }
        }
      }
    }
  }else{stop(paste("No BLAST files found in ",DIR," directory",sep="'"))}
}
combine.fasta=function(AA.sequences,DNA.sequences, AA.recovered.sequences=NULL, SD=2, performMSA=F, Method="ClustalW"){
  
  require(Biostrings)
  require(msa)
  
  if(!is.null(AA.recovered.sequences)){Combined.sequences=c(AA.sequences,AA.recovered.sequences)}else{Combined.sequences=AA.sequences}
  max.AA=0
  for(i in 1:length(Combined.sequences)){L=length(Combined.sequences[[i]]);if(L>max.AA){max.AA=L}}
  write.fasta(Combined.sequences,names=names(Combined.sequences),file.out = paste(Target,'Combined sequences AA.fasta',sep=' '),nbchar = max.AA)
  
  changenames=function(x){unlist(strsplit(x,split=' '))[1]}
  NewNames=as.vector(lapply(names(Combined.sequences),changenames),mode='character')
  Combined.sequences.DNA=DNA.sequences[NewNames]
  max.DNA=0
  for(i in 1:length(Combined.sequences.DNA)){L=length(Combined.sequences.DNA[[i]]);if(L>max.DNA){max.DNA=L}}
  write.fasta(Combined.sequences.DNA,names=names(Combined.sequences.DNA),file.out = paste(Target,'Combined sequences DNA.fasta',sep=' '),nbchar = max.DNA)
  
  if(performMSA){
    AAtoAlign=Biostrings::readAAStringSet(paste(Target,'Combined sequences AA.fasta',sep=' '))
    aln=msa(inputSeqs=AAtoAlign,verbose=T, method=Method)
    dist=as.matrix(seqinr::dist.alignment(msaConvert(aln,"seqinr::alignment")))
    
    mean.dist=mean(rowMeans(dist))
    sd.dist=sd(rowMeans(dist))
    
    outliers.AA=Combined.sequences[rownames(dist[rowMeans(dist)>(mean.dist+SD*sd.dist),])]
    
    if(length(outliers.AA)>0){write.fasta(outliers.AA,names=names(outliers.AA),file.out = paste(Target,'outliers AA.fasta',sep=' '))}
  }
}
