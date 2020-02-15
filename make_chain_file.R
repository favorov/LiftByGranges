#' Create a chain object from a GRanges object genomic coordinates
#' 
#' `GRangesMappingToChainFile` creates a chain object mapped from the genomic coordinates of a GRanges object.
#' 
#' @param input_GRanges A GRanges object that will be converted to a chain file. Intervals must be non-overlapping. Required
#' @param chrom_suffix The suffix to be appended to all chromosome names created in the chain file, if desired.
#' @param verbose Output updates while the function is running. Default FALSE
#' @param alignment The human genome alignment used, either "hg19" or "hg38". Default "hg38"
#' 


GRangesMappingToChainFile<-function(input_GRanges,
                                    chrom_suffix = "",
                                    verbose=FALSE,
                                    alignment="hg38"){
  # first load necessary packages, rtracklayer, plyranges, and TxDb.Hsapiens.UCSC.hg19.knownGene or TxDb.Hsapiens.UCSC.hg38.knownGene
  if("rtracklayer" %in% gsub("package:", "", search())){
    if(verbose==TRUE){print("rtracklayer is loaded correctly")}
  } else {return("please install and/or load rtracklayer")}
  if("plyranges" %in% gsub("package:", "", search())){
    if(verbose==TRUE){print("plyranges is loaded correctly")}
  } else {return("please install and/or load plyranges")}
  
  # determine which genome alignment to use for the chain
  if(!(alignment %in% c("hg19", "hg38"))){
    return("Chromosome alignment must be either hg19 or hg38")
  }
  if(alignment=="hg19"){
    if("TxDb.Hsapiens.UCSC.hg19.knownGene" %in% gsub("package:", "", search())){
      if(verbose==TRUE){print("TxDb.Hsapiens.UCSC.hg19.knownGene is loaded correctly")}
    } else {return("please install and/or load TxDb.Hsapiens.UCSC.hg19.knownGene")}
    seqinft<-as.data.frame(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))
  }
  if(alignment=="hg38"){
    if("TxDb.Hsapiens.UCSC.hg38.knownGene" %in% gsub("package:", "", search())){
      if(verbose==TRUE){print("TxDb.Hsapiens.UCSC.hg38.knownGene is loaded correctly")}
    } else {return("please install and/or load TxDb.Hsapiens.UCSC.hg38.knownGene")}
    seqinft<-as.data.frame(seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))
  }
  
  #create seqinfo object for all chromosomes to get lengths
  seqinf<-as(seqinft[nchar(rownames(seqinft))<6,],"Seqinfo")
  # make certain GRanges chromosome names look like UCSC chromosome names
  if(sum(input_GRanges@seqnames %in% seqinf@seqnames) == 0){
    input_GRanges@seqnames<-paste0("chr", input_GRanges@seqnames) %>% Rle()
    if(sum(input_GRanges@seqnames %in% seqinf@seqnames) == 0){
      return("GRanges chromosomes do not resemble UCSC chromosome names. Suggested format: 'chr#' or just # for chromosome name, e.g. chr1 chr10 chrM")
    }
  }
  # make a chain for each chromosome in the genome
  for(chr in seqinf@seqnames){
    if(verbose==TRUE){print(paste("Chromosome", chr, "starting"))}
    gtf_hold<-input_GRanges %>% filter(seqnames==chr)
    if(length(gtf_hold@ranges)==0){next} # in case of chromosomes without data
    length_chr<-sum(gtf_hold@ranges@width)
    transcripts<- unique(gtf_hold@elementMetadata@listData[["transcript_id"]])
    chrom_chains<-matrix(NA, ncol=3, nrow=(length(gtf_hold)+(2*length(transcripts))))
    row<-1
    start_position=1
    for(t in transcripts){
      gtf_t<-gtf_hold[(elementMetadata(gtf_hold)[,"transcript_id"] == t)]
      length_t<-sum(gtf_t@ranges@width)
      end_position<-start_position+length_t-1
      first_line<-c(paste0("chain 42 ", chr, " ", seqinf@seqlengths[which(seqinf@seqnames==chr)], 
                           " * ", gtf_t@ranges@start[1], " ", 
                           gtf_t@ranges@start[length(gtf_t@ranges@start)]+gtf_t@ranges@width[length(gtf_t@ranges@width)]-1,
                           " ", chr, chrom_suffix, " ", length_chr, " * ", 
                           start_position, " ", end_position," ", t), "", "")
      txpt<-matrix(NA, nrow=length(gtf_t@ranges@start), ncol=3)
      txpt[,1]<-gtf_t@ranges@width
      if(length(gtf_t@ranges@width)>1){ #for cases where there is only 1 exon in a transcript
        txpt[,2]<-c(gtf_t@ranges@start[2:length(gtf_t@ranges@start)]-
                      (gtf_t@ranges@width[1:(length(gtf_t@ranges@width)-1)]+
                       gtf_t@ranges@start[1:(length(gtf_t@ranges@width)-1)])-1, "")
        txpt[,3]<-0
        txpt[nrow(txpt),3]<-""
      }
      chrom_chains[row,]<-first_line
      chrom_chains[(row+1):(row+nrow(txpt)),]<-txpt
      chrom_chains[row+nrow(txpt)+1,]<-c("", "", "")
      start_position<-end_position
      row<-row+nrow(txpt)+2
    }
    # append the chains for the chromosome to the chain file
    write.table(chrom_chains, file=out_chain_name, append=T, sep="\t", row.names = F, col.names = F, quote = F, na="")
    if(verbose==TRUE){print(paste("Chromosome", chr, "complete"))}
  }
  #this is a band-aid fix since import cannot differentiate chains without a truly blank line between chains
  remove_blanks<-readLines(out_chain_name)
  remove_blanks<-gsub("\t\t", "", remove_blanks)
  write(remove_blanks, file=out_chain_name, sep="")
}
