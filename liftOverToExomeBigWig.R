#' Liftover strand-separated BigWigs to an exome BigWig
#' 
#' `liftOverToExomeBigWig` writes a BigWig exome file that is lifted over from two strand-separated BigWig files.
#' 
#' @param input_bw A vector of the names of input BigWig files in the format of c(forward_reads, reverse_reads). If files end in "_forward.bw" and "_reverse.bw", can input the prefix as a single string. Required
#' @param chain The name of the chain file to be used. Required
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
  
}