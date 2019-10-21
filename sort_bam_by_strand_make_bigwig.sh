#!bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -i input_bam -p prefix_output"
   echo -e "\t-i Name of inported BAM file to be turned into strand-separated BAM and BigWig files"
   echo -e "\t-p Prefix of output BAM and BigWig files"
   exit 1 # Exit script after printing help
}

while getopts "i:p:" opt
do
   case "$opt" in
      i ) input_bam="$OPTARG" ;;
      p ) prefix_output="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$prefix_output" ] || [ -z "$input_bam" ]
then
   echo "-i input_bam and -p prefix_output argumetns are required. Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct

samtools view -H $input_bam > header.tmp

samtools view $input_bam | gawk '(and(16,$2))' > forward.tmp

cat header.tmp forward.tmp | samtools view -S -b > "${prefix_output}_forward.bam"

rm forward.tmp

samtools view G3BP1_rep1_sorted.bam | gawk '(! and(16,$2))' > reverse.tmp

cat header.tmp reverse.tmp | samtools view -S -b > "${prefix_output}_reverse.bam"

rm reverse.tmp

rm header.tmp

python bam_to_bigwig_edit.py -i "${prefix_output}_forward.bam" -g hg19 -o "${prefix_output}_forward.bw"

python bam_to_bigwig_edit.py -i "${prefix_output}_reverse.bam" -g hg19 -o "${prefix_output}_reverse.bw"
