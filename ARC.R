ORF.finder=function(sequences,target=NULL,expected.length=NULL,plus=0,minus=0,ANTISENS=F,preferred=NULL,verbose=F){
  T0=Sys.time()
  
  require(seqinr)
  require(dplyr)
  
  if(!is.null(expected.length)){Range=c((expected.length-minus),(expected.length+plus))}else{Range=NULL}
  if(ANTISENS){frames=c(0,1,2,0,1,2)}else{frames=c(0,1,2)}
  
  for (i in 1:length(sequences)){
    
    translation.table=NULL
    for (j in 1:length(frames)){
      if(j<4){S='F'}else{S='R'}
      AA=seqinr::translate(unlist(strsplit(as.character(sequences[[i]]),split='')), frame = frames[j], sens = S)
      
      
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
    
    if(!is.null(preferred)){
      if(preferred=='forward'){if('F'%in%Best.AA$Sens){Best.AA=Best.AA[Best.AA$Sens=='F',]}}
      if(preferred=='reverse'){if('R'%in%Best.AA$Sens){Best.AA=Best.AA[Best.AA$Sens=='R',]}}
    }
    
    if(nrow(Best.AA)==0){writeXStringSet(sequences[i], filepath =  paste(target,"Non ORF DNA.fasta",sep=' '),append = T )}
    
    if(nrow(Best.AA)>=1){
      
      for(n in 1:nrow(Best.AA)){
        AA=seqinr::translate(unlist(strsplit(as.character(sequences[[i]]),split='')), frame = Best.AA$Frame[n], sens = Best.AA$Sens[n])
        if(!is.null(expected.length)){
          if(between(length(AA),min(Range),max(Range))){
            write.fasta(sequences = AA, names = paste(names(sequences)[i],paste('ORF',n,Best.AA$Sens[n],sep='_')), file.out =  paste(target,"ORF AA expected length.fasta",sep=' '), open="a",nbchar = length(AA))
            #only print the reverse complement ORF if no ORF has been found in the sens frame
            if(Best.AA$Sens[n]=='R'){sequences[i]=reverseComplement(sequences[i]);Sens='ReverseComplement'}else{Sens=NULL}
            #write.fasta(sequences = toupper(sequences[[i]]), names = names(sequences)[i], file.out =  paste(paste(target,Sens,sep=''),"ORF DNA expected length.fasta",sep=' '), open="a",nbchar = length(sequences[[i]]))
            writeXStringSet(sequences[i],  filepath =  paste(paste(target,Sens,sep=''),"ORF DNA expected length.fasta",sep=' '), append = T )
          }else{
            write.fasta(sequences = AA, names = paste(names(sequences)[i],paste('ORF',n,Best.AA$Sens[n],sep='_')), file.out =  paste(target,"ORF AA non expected length.fasta",sep=' '), open="a",nbchar = length(AA))
            #only print the reverse complement ORF if no ORF has been found in the sens frame
            if(Best.AA$Sens[n]=='R'){sequences[i]=reverseComplement(sequences[i]);Sens='ReverseComplement'}else{Sens=NULL}
            writeXStringSet(sequences[i], filepath =  paste(paste(target,Sens,sep=''),"ORF DNA non expected length.fasta",sep=' '), append=T)
          }
        }else{
          write.fasta(sequences = AA, names = paste(names(sequences)[i],paste('ORF',n,Best.AA$Sens[n],sep='_')), file.out =  paste(target,"ORF AA.fasta",sep=' '), open="a",nbchar = length(AA))
          #only print the reverse complement ORF if no ORF has been found in the sens frame
          if(Best.AA$Sens[n]=='R'){sequences[i]=reverseComplement(sequences[i]);Sens='ReverseComplement'}else{Sens=NULL}
          writeXStringSet(sequences[i], filepath =  paste(paste(target,Sens,sep=''),"ORF DNA.fasta",sep=' '),append=T)
        }
      }
    }
    
    if(verbose==T){print( paste('Processing: ',round((i/length(sequences))*100,4),'%',sep='') )}
  }
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
