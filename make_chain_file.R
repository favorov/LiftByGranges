#' Create exome chain file from genome annotation file
#' 
#' `GRangesMappingToChainFile` writes a chain file of the exome mapped from a genome annotation file.
#' 
#' @param input_GRange The name of the gtf/gff file that will be converted to an exome chain file. Required
#' @param out_chain_name The name of the chain file to be created. Required
#' @param chrom_suffix The suffix to be appended to all chromosome names created in the chain file, if desired.
#' @param verbose Output updates while the function is running. Default FALSE
#' @param transcript_list a vector of transcript names that represent the most expressed isoform of their respective genes and correspond to gtf annotation names. Required
#' @param alignment The human genome alignment used, either "hg19" or "hg38". Default "hg19"
#' 
#' @examples 
#' \dontrun{
#' GRangesMappingToChainFile("hg19_annotations_ensembl.gtf",
#'   "hg19_annotations_ensembl_exome.chain",
#'   c("ENST00000600805", "ENST00000595082", "ENST00000597230", "ENST00000595005")
#'   chrom_suffix="exome_muscle",
#'   verbose=T)
#' }
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
  
  # get alignment info, make sure it is compatible
  if(!(alignment %in% c("hg19", "hg38"))){
    return("valid genome alignments are hg19 and hg38")
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
  # write chain file chromosome-by-chromosome
  if(sum(input_GRanges@seqnames %in% seqinf@seqnames) == 0){
    input_GRanges@seqnames<-paste0("chr", input_GRanges@seqnames) %>% Rle()
    if(sum(input_GRanges@seqnames %in% seqinf@seqnames) == 0){
      return("GTF chromosomes do not resemble UCSC chromosome names. Suggested format: 'chr#' or just # for chromosome name, e.g. chr1 chr10 chrM")
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
  chrom_chains<-
  test<-as.character(chr)
    # append the chains for the chromosome to the chain file
    #write.table(chrom_chains, file=out_chain_name, append=T, sep="\t", row.names = F, col.names = F, quote = F, na="")
    #if(verbose==TRUE){print(paste("Chromosome", chr, "complete"))}

  #this is a band-aid fix since import cannot differentiate chains without a truly blank line between chains
  #remove_blanks<-readLines(out_chain_name)
  #remove_blanks<-gsub("\t\t", "", remove_blanks)
  #write(remove_blanks, file=out_chain_name, sep="")
}
