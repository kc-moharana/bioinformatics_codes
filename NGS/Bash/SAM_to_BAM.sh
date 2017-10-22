
## Shell script to convert SAM to sorted and indexed BAM file
## Change threads manually 
## Shell script to convert SAM to sorted and indexed BAM file

if [ $# -lt 2 ]
then
	echo "Shell script to convert SAM to sorted and indexed BAM file"
	echo "Usage: bash $0 SAMFILE [CPU]"
  exit 1
fi

SAM=$1
threads=$2 

echo "samtools view -bh -@ $threads -o $SAM\.bam $SAM"
samtools view -bh -@ $threads -o $SAM\.bam $SAM
#rm $SAM
samtools sort -@ $threads -o $SAM\.sort\.bam $SAM\.bam
rm $SAM\.bam

echo "samtools index $SAM\.sort\.bam "
samtools index $SAM\.sort\.bam 

echo "[[MESSAGE]] :: Delete original SAM file : $SAM"



