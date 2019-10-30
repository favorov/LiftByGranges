#' Make a FASTA list of transcript sequences from a genome FASTA file, genome annotation file, and list of transcript names
#' 
#' `ExtractTranscriptomeSequence` writes a FASTA file of transcript sequences from a list of transcripts.
#' 
#' @param transcript_list a vector of transcript names that represent the most expressed isoform of their respective genes and correspond to GTF annotation names. Required
#' @param ref_genome The name of the reference genome FASTA from which exome sequences will be derived. Required
#' @param transcript_gtf The name of the exome GTF annotation file. Coordinates must match the file input for the ref_genome parameter. Required if there is no genome_gtf
#' @param genome_gtf The name of the GTF/GFF file that contains all exome annotations. Coordinates must match the file input for the ref_genome parameter. Required if there is no transcript_gtf
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

ExtractTranscriptomeSequence<-function(transcript_list, ref_genome, transcript_gtf,
                                       genome_gtf, exome_prefix="exome"){
  # load necessary packages plyranges and rtracklayer
  if("rtracklayer" %in% gsub("package:", "", search()) == F){print("please install and/or load rtracklayer")}
  if("plyranges" %in% gsub("package:", "", search()) == F){print("please install and/or load plyranges")}
  
  # make sure there is a gtf file available; if there is, load in the information; if not, make one
  if(missing(genome_gtf)){
    if(missing(transcript_gtf)){
      return("Please provide either a genome or a transcriptome GTF file via the genome_gtf or transcript_gtf parameters.")
    } else{ input_gtf <- transcript_gtf}
  } else{
    gtf<-import(genome_gtf)
    gtf<-gtf %>% filter(type == "exon") #since only care about exome
    gtf_transcripts<-gtf[(elementMetadata(gtf)[,"transcript_id"] %in% transcript_list)]
    input_gtf<-paste0(exome_prefix, ".gtf")
    write_gff(gtf_transcripts, input_gtf)
  }
  
  # call bedtools on local system, get FASTA sequences
  system(paste0("bedtools getfasta -s -fi ", ref_genome, " -bed ", input_gtf, " -fo ", exome_prefix, ".fa"))

}
