#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my $usage = <<usage;
        Usage:perl $0 -ref -bam -dup -chr -pos |-list
	    -ref    reference sequence
	    -bam    alignment file
	    -dup    whether or not stat PCR duplicate reads,default delete dup
	    -chr    chromosome
	    -pos    position
	    -check  base, output readsID with this base
	    -list   file of positions list, parallel with -chr -pos combined
        Example:
	    perl $0 -ref -bam -chr -pos
	    perl $0 -ref -bam -list
       note:

	list format: chr pos(1 base) 
	samtools should be installed
usage

my ( $bam, $chr, $pos, $check, $list, $rev, $dup, $out, $fa, $ref,$aim );
GetOptions(
    "bam:s"    => \$bam,
    "chr:s"    => \$chr,
    "pos:s"    => \$pos,
    "check:s"   => \$check,
    "list:s"   => \$list,
    "reverse!" => \$rev,
    "dup!"     => \$dup,
    "obam:s"   => \$out,
    "ref:s"    => \$ref,
    "fa:s"    => \$fa,
    "aim:s"    => \$aim,
);
die $usage
  unless ( defined $bam
    && ( defined $list || ( defined $chr && defined $pos ) ) );
if(-d "/ifs7"){
$ref||="/ifs7/B2C_SGD/PROJECT/PP12_Project/analysis_pipeline/HPC_chip/db/aln_db/hg19/hg19_chM_male_mask.fa";
}elsif(-d "/jdfstj2"){
$ref="/jdfstj2/B2C_RD_P2/backup/share/udata/wangyaoshen/src/HPC_chip/db/aln_db/hg19/hg19_chM_male_mask.fa";
}

#my$samtools="/ifs9/BC_B2C_01A/B2C_SGD/Newborn/analysis_pipeline/HPC_chip/tools/samtools-0.1.19";
#$ref||="/hwfssz1/ST_MCHRI/CLINIC/SOFTWARES/analysis_pipeline/HPC_chip/db/aln_db/hg19/hg19_chM_male_mask.fa";
open F,">$fa" or die "can't output fa\n";
my$samtools="samtools";
my $sample=basename( $bam);
if ( defined $list ) {
    open L, "zcat -f $list|" or die $!;
    while (<L>) {
        chomp;
        my ( $cc, $pp ) = split;
        parse( $bam, $cc, $pp );
    }
}
else {
    parse( $bam, $chr, $pos );
}

sub parse {
    my ( $bam, $chr, $pos ) = @_;
    $chr =~ s/chr//ig;
    $chr = "chr$chr";

    my ( $pos1, $pos2, $type, $base, %fa, %stat );

    if ( $pos =~ /^(\d+)([ATCG])?$/ ) {
        $pos1 = $pos2 = $1;
        $type = "mis";
        $base = $2;
    }
    elsif ( $pos =~ /^(\d+)in(\d+)([ATCG]+)$/ ) {
        $pos1 = $1;
        $pos2 = $2;
        $type = "ins";
        $base = $3;
    }
    elsif ( $pos =~ /^(\d+)to(\d+)$/ ) {
        $pos1 = $1;
        $pos2 = $2;
        $type = "del";

    }
    $base||="";
    if ( defined $dup ) {
        open B, "$samtools view $bam $chr:$pos1-$pos2|" or die $!;
    }
    else { open B, "$samtools view -F 0x400 $bam $chr:$pos1-$pos2|" or die $!; }
    my$total_depth;
    while (<B>) {
        chomp;
        my @linetmp = split(/\t/);
        my ( $flag, $rname, $pos, $mq, $cigar, $seq, @extras ) =
          @linetmp[ 1 .. 5, 9, 11 .. $#linetmp ];
        my @num = split /\D/,  $cigar;    #105##29##
        my @mat = split /\d+/, $cigar;    #M##S##
        shift @mat;

        # print scalar @num,"\n",scalar @mat,"\n";
        my ( $start, $end );
        for my $i ( 0 .. $#num ) {
            if ( $mat[$i] =~ /M/ ) {
                #	    print $i."\n";
                if ( $i == 0 ) {
                    $start = $pos;
                    $end += ( $start + $num[$i] - 1 );
                }
                else {
                    #if(!defined $end){print $linetmp[0]."\n";exit;}
                    $start = $end + 1;
                    $end += $num[$i];
                }
                #	print join "\t",$start,$end,$pos1,"\n";
                if (   $type eq "mis"
                    && $rname eq $chr
                    && ( $pos1 >= $start && $pos1 <= $end ) )
                {    ##only compare  mismatch type,and in the range
                        #    my @refseq = split /\n/,
                        # `samtools faidx $ref $rname:$pos1-$pos1`;
                        #shift @refseq;
                        #my $refseq = join "", @refseq; no use
                    my $dis = $pos1 - $start;
                    my $base2 = substr( $seq, $dis, 1 ); ##����bam����е����ж��Ѿ��Ͳο�����һ���ˣ���Ҫ���컥���ľͷ��򻥲�������
                    if ( defined $check && $base2 eq $check ) {
                        #$fa{ $linetmp[0] } = 0;    ##mark reads ID
                        print F "$linetmp[0]\n";
                    }


		    $total_depth++;
			#if($base2 eq $aim){
			#	print $_."\n";				
			#}
		    if($flag & 0x10 ){
			$stat{$base2}{"+"}++;
			}else{
			$stat{$base2}{"-"}++;
 			}##strand record

                }
                $seq = substr( $seq, $num[$i] );    ##sequence delete cycly
            }
            elsif ( $mat[$i] =~ /D/ ) {             ##deletion
                if ( $i == 0 ) {
                    $start = $pos;
                    $end += ( $start + $num[$i] - 1 );
                }
                else {
                    $start = $end + 1;
                    $end += $num[$i];
                }
                ##deletion compare
                if ( $start == $pos1 && $end == $pos2 ) {
                    $fa{ $linetmp[0] } = 0;
                    ##$seq=substr($seq,$num[$i]);	##no sequence in the reads
                }
            }
            elsif ( $mat[$i] =~ /I/ ) {    ##insertion
                                           #133M3I10M5S
                my $base2 = substr( $seq, 0, $num[$i] );
                if ( $end == $pos1 && $base2 eq $base ) {
                    $fa{ $linetmp[0] } = 0;
                }
                $seq = substr( $seq, $num[$i] );    ##sequence delete cycly
            }
            elsif ( $mat[$i] =~ /S/ ) {             ##soft clip 144M7S  7S144M
                if ( $i == 0 ) {
                    $start = $pos - 1;
                    $end   = $pos - 1;
                }
                $seq = substr( $seq, $num[$i] );    ##sequence delete cycly
            }
            elsif ( $mat[$i] =~ /H/ ) {
                if ( $i == 0 ) {
                    $start = $pos - 1;
                    $end   = $pos - 1;
                }
            }

        }
    }
    close B;
    if ( $type eq "mis" ) {
	my$rr=`$samtools faidx $ref $chr:$pos1-$pos1|awk 'NR==2'`;chomp $rr;
        print "$sample\t$chr\t$pos1\t$rr";
	#$total_depth
        for my $i ( "A", "T", "C", "G","N" ) {
            if ( exists $stat{$i} ) {
		my$tg;
		my$tt;
		if(exists $stat{$i}{"+"}){$tg.=$stat{$i}{"+"}."+,";$tt+=$stat{$i}{"+"};}else{$tg.="0+,";}
		if(exists $stat{$i}{"-"}){$tg.=$stat{$i}{"-"}."-";$tt+=$stat{$i}{"-"};}else{$tg.="0-";}
		my$rt=sprintf("%.2f",$tt/$total_depth*100);
		$rt="$rt%";
                print "\t$i:" . $tt."($rt,$tg)";
            }
            else {
                print "\t$i:0";
            }
        }
        print "\n";
    }
}
close F;
## Please see file perltidy.ERR
