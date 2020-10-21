#am=../bam/19B2767722.ref_CYP21A2.final.bam
#os=32007584
#am=19B2767722
#ase=A
bam=$1
sam=$2
pos=$3
base=$4
alt=$5
if [[ $# -lt 5 ]];
then
	echo sh $0 bam samid pos base ref_or_alt_tag
exit
fi

tag=$sam.$pos.$alt$base
perl /zfsyt1/B2C_RD_P2/PMO/P18Z12204N0372/CYP21A2/shell/support_reads.d.pl -ref $hg19 -bam $bam -chr chr6 -pos $pos -check $base -fa $tag.fa

IFS=$'\n' read -d '' -r -a lines < $tag.fa
bar=$(printf "|%s" "${lines[@]}")
bar=${bar:1}
samtools view -H $bam >hd.sam
st=$(($pos-2000))
ed=$(($pos+2000))
samtools view $bam chr6:$st-$ed|egrep -w "$bar" >$tag.sam
cat hd.sam $tag.sam|samtools view -bho $tag.bam -
samtools index $tag.bam
rm $tag.fa hd.sam $tag.sam

