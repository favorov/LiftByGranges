#' Create exome chain file from genome annotation file
#' 
#' `GRangesMappingToChainFile` writes a chain file of the exome mapped from a genome annotation file.
#' 
#' @param input_gtf The name of the gtf/gff file that will be converted to an exome chain file
#' @param out_chain_name The name of the chain file to be created
#' @param chrom_suffix The suffix to be appended to all chromosome names created in the chain file
#' 

GRangesMappingToChainFile<-function(input_gtf, out_chain_name, chrom_suffix = "exome"){
  # first install necessary packages, rtracklayer and plyranges
  if(require("rtracklayer")){
    print("rtracklayer is loaded correctly")
  } else {
    print("trying to install rtracklayer")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
      BiocManager::install("rtracklayer")
    if(require(rtracklayer)){
      print("rtracklayer installed and loaded")
    } else {
      stop("could not install rtracklayer")
    }
  }
  if(require("plyranges")){
    print("plyranges is loaded correctly")
  } else {
    print("trying to install plyranges")
    install.packages("plyranges")
    if(require(plyranges)){
      print("plyranges installed and loaded")
    } else {
      stop("could not install plyranges")
    }
  }
  # delete old chain file if file with this name already exists
  system(paste0("test -e ", out_chain_name, " && rm ", out_chain_name))
  gtf<-import(input_gtf)
  gtf<-gtf %>% filter(type == "exon") #since only care about post-processed transcriptome
  lapply(as.character(gtf@seqnames@values), function(chr){
    for(str in c("-", "+")){
      gtf_hold<-gtf %>% filter(strand==str, seqnames==chr)
      
      if(length(gtf_hold@ranges)==0){next} # in case of chromosomes without data on one strand
      length_chr<-sum(gtf_hold@ranges@width)
      if(str=="+"){sign<-"plus"}else{sign<-"minus"}
      transcripts<- unique(gtf_hold@elementMetadata@listData[["transcript_id"]])
      
      ####
      arbitrary_length<-max(gtf_hold@ranges@start)+ max(gtf_hold@ranges@width) #### NEED TO GET ACTUAL CHROMOSOME LENGTHS
      ####
      
      # need to turn this into a table that then gets written all at once
      # instead of thousands of sequencial appendations
      
      n=1
      for(t in transcripts){
        gtf_t<-gtf_hold %>% filter(transcript_id == t)
        length_t<-sum(gtf_t@ranges@width)
        n2<-n+length_t-1
        first_line<-paste0("chain 42 ", chr, "_", sign, "_", chrom_suffix, " ", length_chr, " ", str,
                           " ", n, " ", n2, " ", chr, " ", arbitrary_length, " ",  gtf_t@ranges@start[1], " ",
                           gtf_t@ranges@start[length(gtf_t@ranges@start)]+gtf_t@ranges@width[length(gtf_t@ranges@width)]-1,
                           " ", t)
        
        txpt<-matrix(NA, nrow=length(gtf_t@ranges@start), ncol=3)
        txpt[,1]<-gtf_t@ranges@width
        if(length(gtf_t@ranges@width)>1){ #for cases where there is only 1 exon in a transcript
          txpt[,2]<-0
          txpt[nrow(txpt),2]<-""
          txpt[,3]<-c(gtf_t@ranges@start[2:length(gtf_t@ranges@start)]-
                        (gtf_t@ranges@width[1:(length(gtf_t@ranges@width)-1)]+
                           gtf_t@ranges@start[1:(length(gtf_t@ranges@width)-1)])-1, "")
        }
        
        write(paste0("\n", first_line), file=out_chain_name, append=T)
        write.table(txpt, file=out_chain_name, append=T, sep="\t", row.names = F, col.names = F)
        
        n<-n2
      }
    }
  }) %>% invisible()
}




