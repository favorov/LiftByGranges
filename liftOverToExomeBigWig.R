#' Liftover strand-separated BigWigs to an exome BigWig
#' 
#' `liftOverToExomeBigWig` writes a BigWig exome file that is lifted over from two strand-separated BigWig files.
#' 
#' @param input_bw A vector of the names of input BigWig files in the format of c(forward_reads, reverse_reads). If file names end in '_forward.bw' and '_reverse.bw', it can be a single string corresponding to the prefix of the files. Required
#' @param chain The name of the chain file to be used. Formet should be like chain files derived from `GRangesMappingToChainFile`. Required
#' @param output_bw The name of the exome BigWig file. Required
#' @param write_chr Whether to output a table of mapped chromosome names and lengths. Default TRUE
#' @param out_chr Name of the chromosome names and lengths table file. Required if `write_chr` is TRUE
#' 
#' @examples 
#' \dontrun{
#' liftOverToExomeBigWig("CLIP_rep1",
#' "hg19_exome.chain",
#' "CLIP_exome.bw",
#' "chain.chr.size")
#' 
#' #' liftOverToExomeBigWig(c("CLIP_rep1_F.bw", "CLIP_rep1_R.bw"),
#' "hg19_exome.chain",
#' "CLIP_exome.bw",
#' write_chr=F)
#' }
#' 

liftOverToExomeBigWig<-function(input_bw, chain, output_bw, write_chr = TRUE, out_chr){
  # first load necessary packages, rtracklayer, plyranges, and dplyr
  if("rtracklayer" %in% gsub("package:", "", search()) == F){print("please install and/or load rtracklayer")}
  if("plyranges" %in% gsub("package:", "", search()) == F){print("please install and/or load plyranges")}
  if("dplyr" %in% gsub("package:", "", search()) == F){print("please install and/or load dplyr")}

  if(length(input_bw)==2){
    f_bw<-input_bw[1]
    r_bw<-input_bw[2]
  } else if(length(input_bw)==1){
    f_bw<-paste0(input_bw, "_forward.bw")
    r_bw<-paste0(input_bw, "_reverse.bw")
  } else{print("input_bw must be in in the format of c(forward_reads, reverse_reads). If file names end in '_forward.bw' and '_reverse.bw', input_bw can be a single string corresponding to the prefix of the files")}
  
  #since reads are rev. complements, assign reads to the opposite strand
  plus<-import(f_bw)
  strand(plus)<-"+"
  minus<-import(r_bw)
  strand(minus)<-"-"
  
  chain_object<-import(chain, exclude="junk")
  
  liftOver_plus<-liftOver(plus, chain_object)
  # remove any blank intervals (shouldn't be any >1)
  liftOver_plus<-liftOver_plus[(elementNROWS(range(liftOver_plus))==1L)] %>% unlist()
  liftOver_plus<-liftOver_plus[lapply(as.character(liftOver_plus@seqnames),
                                                          function(x){grepl("plus",x) == T}) %>% unlist()]
  liftOver_minus<-liftOver(minus, chain_object)
  # remove any blank intervals (shouldn't be any >1)
  liftOver_minus<-liftOver_minus[(elementNROWS(range(liftOver_minus))==1L)] %>% unlist()
  liftOver_minus<-liftOver_minus[lapply(as.character(liftOver_minus@seqnames),
                                                            function(x){grepl("minus",x) == T}) %>% unlist()]
  
  liftOver<-bind_ranges(liftOver_minus, liftOver_plus)
  
  #must re-introduce chromosome lengths into seqinfo
  lines <- c()
  for ( chr in liftOver@seqinfo@seqnames ) {
    con <- file(chain, "r")
    while(TRUE) {
      line = readLines(con, 1)
      if(length(line) == 0) break
      else if(grepl(chr, line)){
        lines <- c(lines, line)
        break
      }
    }
    close(con)
  }
  seq.info<-do.call(rbind.data.frame, lapply(lines, function(x){strsplit(x, " ") %>% unlist()}))[,8:9]
  colnames(seq.info)<-c("chr", "length")
  if(write_chr==T){
    #save seq.length file
    write.table(seq.info, out_chr, sep="\t", col.names = F, row.names = F, quote = F)
  }
  liftOver@seqinfo@seqlengths<-seq.info[,2] %>%
    as.character() %>% as.integer()
  liftOver@elementMetadata$score<-1:length(liftOver) #arbitrary score
  liftOver<-liftOver[countOverlaps(liftOver,liftOver) <= 1L] # remove duplicates
  
  rtracklayer::export(liftOver, con= output_bw)
}


