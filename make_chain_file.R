#' GRangesMappingToChainViaFile
#' 
#' GRangesMappingToChainViaFile creates a \code{Chain} object based on intervals from a \code{GRanges} object. The mapping by this \code{Chain} will collapse chromosomes of the annotation into pseudogenome that are combined from tiled intervals of the \code{GRanges} object. 
#' 
#' @param input_GRanges A GRanges file with non-overlapping intervals that will be converted to a chain file. Required.
#' @param out_chain_name The name of a chain file to be written in the local directory. Optional.
#' @param chrom_suffix The suffix to be appended to all chromosome names created in the chain file, if desired.
#' @param verbose Output updates while the function is running. Default FALSE
#' @param alignment The human genome alignment used, either "hg19" or "hg38". Default "hg38"
#' @return a Chain object that maps all the chromosomes according to GRanges

GRangesMappingToChainViaFile<-function(input_GRanges,
                                    chrom_suffix = "",
                                    out_chain_name = "",
                                    verbose=FALSE,
                                    alignment="hg38"){
  # first load necessary packages, rtracklayer, plyranges, and TxDb.Hsapiens.UCSC.hg19.knownGene or TxDb.Hsapiens.UCSC.hg38.knownGene
  if(!("rtracklayer" %in% gsub("package:", "", search()))){stop("please install and/or load rtracklayer")}
  if(!("plyranges" %in% gsub("package:", "", search()))){stop("please install and/or load plyranges")}
  
  # get alignment info, make sure it is compatible
  if(!(alignment %in% c("hg19", "hg38"))){
    stop("valid genome alignments are hg19 and hg38")
  }
  if(alignment=="hg19"){
    if("TxDb.Hsapiens.UCSC.hg19.knownGene" %in% gsub("package:", "", search())){
      if(verbose==TRUE){print("TxDb.Hsapiens.UCSC.hg19.knownGene is loaded correctly")}
    } else {stop("please install and/or load TxDb.Hsapiens.UCSC.hg19.knownGene")}
    seqinft<-as.data.frame(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))
  }
  if(alignment=="hg38"){
    if("TxDb.Hsapiens.UCSC.hg38.knownGene" %in% gsub("package:", "", search())){
      if(verbose==TRUE){print("TxDb.Hsapiens.UCSC.hg38.knownGene is loaded correctly")}
    } else {stop("please install and/or load TxDb.Hsapiens.UCSC.hg38.knownGene")}
    seqinft<-as.data.frame(seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))
  }
  
  #confirm GRanges doesn't have any overlapping intervals
  if(! identical(input_GRanges,reduce(input_GRanges))){
    stop("GRanges object includes overlapping intervals. This will cause errors when trying to use the chain object.")
  }
  #create seqinfo object for all chromosomes to get lengths
  seqinf<-as(seqinft[nchar(rownames(seqinft))<6,],"Seqinfo")
  # write chain file chromosome-by-chromosome
  if(sum(input_GRanges@seqnames %in% seqinf@seqnames) == 0){
    input_GRanges@seqnames<-paste0("chr", input_GRanges@seqnames) %>% Rle()
    if(sum(input_GRanges@seqnames %in% seqinf@seqnames) == 0){
      stop("GTF chromosomes do not resemble UCSC chromosome names. Suggested format: 'chr#' or just # for chromosome name, e.g. chr1 chr10 chrM")
    }
  }
  
  row<-1
  interval<-1
  chrom_chains<-matrix(NA, ncol=3, nrow=(length(input_GRanges)))
  
  # make a chain for each chromosome in the genome
  for(chr in seqinf@seqnames){
    if(verbose==TRUE){print(paste("Chromosome", chr, "starting"))}
    gtf_hold<-input_GRanges %>% filter(seqnames==chr)
    if(length(gtf_hold@ranges)==0){next} # in case of chromosomes without data
    length_chr<-sum(gtf_hold@ranges@width)
    first_line<-c(paste0("chain 42 ", chr, " ", seqinf@seqlengths[which(seqinf@seqnames==chr)], 
                         " * ", gtf_hold@ranges@start[1], " ", 
                         gtf_hold@ranges@start[length(gtf_hold@ranges@start)]+gtf_hold@ranges@width[length(gtf_hold@ranges@width)]-1,
                         " ", chr, chrom_suffix, " ", length_chr, " * 1 ", 
                         length_chr," ", interval), "", "")
    interval<-interval + 1
    chrom<-matrix(NA, nrow=length(gtf_hold@ranges@start), ncol=3)
    chrom[,1]<-gtf_hold@ranges@width
    if(length(gtf_hold@ranges@width)>1){ #for cases where there is only 1 interval in a chromosome
      chrom[,2]<-c(gtf_hold@ranges@start[2:length(gtf_hold@ranges@start)]-
                  (gtf_hold@ranges@width[1:(length(gtf_hold@ranges@width)-1)]+
                  gtf_hold@ranges@start[1:(length(gtf_hold@ranges@width)-1)])-1, "")
      chrom[,3]<-0
      chrom[nrow(chrom),3]<-""
      } 
    chrom_chains[row,]<-first_line
    chrom_chains[(row+1):(row+nrow(chrom)),]<-chrom
    chrom_chains[row+nrow(chrom)+1,]<-c("", "", "")
    row<-row+nrow(chrom)+2
  }
  chrom_chains<-na.omit(chrom_chains)
  if(verbose == TRUE){print("Creating chain object")}
  format_chrom_chains<-apply(chrom_chains, 1, paste, collapse = "\t")
  format_chrom_chains<-gsub("\t\t", "", format_chrom_chains)
  tmp<-""
  if(out_chain_name == ""){
    out_chain_name<-tempfile(pattern = "", fileext = ".chain")
    tmp<-out_chain_name
  } else{if(verbose == TRUE){print("Saving chain object")}}
  writeLines(format_chrom_chains, con=out_chain_name)
  # have to write intermediate file to inport chain object
  chain<-import.chain(out_chain_name)
  unlink(tmp)
  return(chain)
}
