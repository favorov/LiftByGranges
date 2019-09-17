#' Create exome chain file from genome annotation file
#' 
#' `GRangesMappingToChainFile` writes a chain file of the exome mapped from a genome annotation file.
#' 
#' @param input_gtf The name of the gtf/gff file that will be converted to an exome chain file
#' @param out_chain_name The name of the chain file to be created
#' @param chrom_suffix The suffix to be appended to all chromosome names created in the chain file
#' @param verbose Output updates while the function is running. Default FALSE
#' @param transcript_list a vector of transcript names that represent the most expressed isoform of their respective genes and correspond to gtf annotation names
#' 
#' @examples 
#' \dontrun{
#' GRangesMappingToChainFile("hg19_annotations_ensembl.gtf",
#'   "hg19_annotations_ensembl_exome.chain",
#'   c("ENST00000600805", "ENST00000595082", "ENST00000597230", "ENST00000595005")
#'   chrom_suffix="exome_muscle",
#'   verbose=T,)
#' }
#' 

GRangesMappingToChainFile<-function(input_gtf, out_chain_name, transcript_list, chrom_suffix = "exome", verbose=FALSE){
  # first load necessary packages, rtracklayer, plyranges, and TxDb.Hsapiens.UCSC.hg19.knownGene
  if(require("rtracklayer")){
    if(verbose==TRUE){print("rtracklayer is loaded correctly")}
  } else {print("please install rtracklayer")}
  if(require("plyranges")){
    if(verbose==TRUE){print("plyranges is loaded correctly")}
  } else {print("please install plyranges")}
  if(require("TxDb.Hsapiens.UCSC.hg19.knownGene")){
    if(verbose==TRUE){print("TxDb.Hsapiens.UCSC.hg19.knownGene is loaded correctly")}
  } else {print("please install TxDb.Hsapiens.UCSC.hg19.knownGene")}
  
  #create seqinfo object for all chromosomes to get lengths
  seqinft19<-as.data.frame(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))
  seqinf19<-as(seqinft19[nchar(rownames(seqinft19))<6,],"Seqinfo")
  # delete old chain file if file with this name already exists
  system(paste0("test -e ", out_chain_name, " && rm ", out_chain_name))
  gtf<-import(input_gtf)
  gtf<-gtf %>% filter(type == "exon") #since only care about exome
  gtf_transcripts<-gtf[(elementMetadata(gtf)[,"transcript_id"] %in% transcript_list)]
  if(verbose==TRUE){print("annotation data finished loading")}
  for(chr in seqinf19@seqnames){
    for(str in c("-", "+")){
      if(verbose==TRUE){print(paste("Chromosome", chr, "strand", str, "starting"))}
      gtf_hold<-gtf_transcripts %>% filter(strand==str, seqnames==chr)
      if(length(gtf_hold@ranges)==0){next} # in case of chromosomes without data on one strand
      length_chr<-sum(gtf_hold@ranges@width)
      if(str=="+"){sign<-"plus"}else{sign<-"minus"}
      transcripts<- unique(gtf_hold@elementMetadata@listData[["transcript_id"]])
      chrom_chains<-matrix(NA, ncol=3, nrow=(length(gtf_hold)+(2*length(transcripts))))
      row<-1
      start_position=1
      # make a chain for each transcript in the chromosome
      for(t in transcripts){
        gtf_t<-gtf_hold[(elementMetadata(gtf_hold)[,"transcript_id"] == t)]
        length_t<-sum(gtf_t@ranges@width)
        end_position<-start_position+length_t-1
        first_line<-c(paste0("chain 42 ", chr, " ", seqinf19@seqlengths[which(seqinf19@seqnames==chr)], 
                             " ", str, " ", gtf_t@ranges@start[1], " ", 
                             gtf_t@ranges@start[length(gtf_t@ranges@start)]+gtf_t@ranges@width[length(gtf_t@ranges@width)]-1,
                             " ", chr, "_", sign, "_", chrom_suffix, " ", length_chr, " ", str, " ", 
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
      if(verbose==TRUE){print(paste("Chromosome", chr, "strand", str, "complete"))}
    }
  }
  #this is a band-aid fix since import cannot differentiate chains without a truly blank line between chains
  remove_blanks<-readLines(out_chain_name)
  remove_blanks<-gsub("\t\t", "", remove_blanks)
  write(remove_blanks, file=out_chain_name, sep="")
}
