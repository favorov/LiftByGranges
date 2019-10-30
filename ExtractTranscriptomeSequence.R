#' Make a FASTA list of transcript sequences from a genome FASTA file, genome annotation file, and list of transcript names
#' 
#' `ExtractTranscriptomeSequence` writes a FASTA file of transcript sequences from a list of transcripts.
#' 
#' @param transcript_list a vector of transcript names that represent the most expressed isoform of their respective genes and correspond to GTF annotation names. Required
#' @param ref_genome The name of the reference genome FASTA from which exome sequences will be derived. Required
#' @param genome_gtf The name of the GTF/GFF file that contains all exome annotations. Coordinates must match the file input for the ref_genome parameter. Required
#' @param exome_prefix The prefix for all exome output files, including the FASTA. Default "exome"
#' @param bedtools_path Necessary for individuals running R in a GUI on some versions of OSX
#' 
#' @examples 
#' \dontrun{
#' ExtractTranscriptomeSequence(c("ENST00000600805", "ENST00000595082", "ENST00000597230", "ENST00000595005"),
#'   ref_genome = "hg19.fa",
#'   genome_gtf = "hg19_annotations_ensembl.gtf")
#' }
#'

ExtractTranscriptomeSequence<-function(transcript_list, ref_genome,
                                       genome_gtf, exome_prefix="exome"){
  # load necessary packages plyranges and rtracklayer
  if("rtracklayer" %in% gsub("package:", "", search()) == F){print("please install and/or load rtracklayer")}
  if("plyranges" %in% gsub("package:", "", search()) == F){print("please install and/or load plyranges")}
  
  # make sure there is a gtf file available; if there is, load in the information; if not, make one
  if(missing(genome_gtf)){
      return("Please provide a genome GTF file via the genome_gtf parameter.")
  }else{
      gtf<-import(genome_gtf)
      gtf<-gtf %>% filter(type == "exon") #since only care about exome
      gtf_transcripts<-gtf[(elementMetadata(gtf)[,"transcript_id"] %in% transcript_list)]
      input_bed<-paste0(exome_prefix, ".bed")
      df <- data.frame(seqnames=seqnames(gtf_transcripts),
                       starts=start(gtf_transcripts)-1,
                       ends=end(gtf_transcripts),
                       names=elementMetadata(gtf_transcripts)[,"transcript_id"],
                       strands=strand(gtf_transcripts))
      write.table(df, input_bed, sep="\t", quote = F, row.names = F, col.names = F)
  }
  
  # call bedtools on local system, get FASTA sequences
  system(paste0("bedtools getfasta -s -fi ", ref_genome, " -bed ", input_bed, " -name -tab -fo ", exome_prefix, ".tmp"))

  # BUT FASTA are fragmented; must compile sequences into whole transcripts to input into folding algorithm
  edit_tbl<-read.table(paste0(exome_prefix, ".tmp"), stringsAsFactors = F)
  system(paste0("rm ", exome_prefix, ".tmp"))
  edit_tbl$V1<-substr(edit_tbl$V1, 1, nchar(edit_tbl$V1)-2)
  edit_tbl$strand<-df$strand[match(df$names, edit_tbl$V1)] #get strand info
  
  # write.fasta function is simplified from seqinr write.fasta
  write.fasta <- function(sequences, names, file.out){
    # Open output file:
    outfile <- file(description = file.out, open = "w")
    # Function to write one sequence in output file:
    write.oneseq <- function(sequence, name){
      writeLines(paste(">", name, sep = ""), outfile)
      writeLines(sequence, outfile)
    }
    # Write all sequences in output file:
      n.seq <- length(sequences)
      sapply(seq_len(n.seq), function(x) write.oneseq(sequence = as.character(sequences[[x]]), name = names[x]))
      close(outfile)
  }
  
  hold_matrix<-matrix(NA, ncol=2, nrow=length(unique(edit_tbl$V1)))
  hold_matrix[,1]<-unique(edit_tbl$V1)
  for(txpt in unique(edit_tbl$V1)){
    filter(edit_tbl, V1 == txpt) -> txpt_tbl
    if(unique(txpt_tbl$strand)== "-"){
      seq<-paste(rev(txpt_tbl$V2), collapse="", sep="")
    }else{
      seq<-paste(txpt_tbl$V2, collapse="", sep="")
    }
    hold_matrix[which(hold_matrix[,1]==txpt),2]<-seq
  }
  
  write.fasta(hold_matrix[,2], hold_matrix[,1], paste0(exome_prefix, ".fa"))
}
