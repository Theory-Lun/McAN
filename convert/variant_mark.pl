#!/usr/bin/perl -w
#Author: tiandm
#Date: 2022-01-21
#Description: 输入多序列比对文件，指定参考序列ID，鉴定序列中变异位点数据

my $version=8.00;

use strict;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);

my %opts;
GetOptions(\%opts,"fas=s","ref=s","o=s","h");
if (!(defined $opts{fas} and defined $opts{ref} and defined $opts{o}) || defined $opts{h}) { #necessary arguments
&usage;
}

#read in aligne data
my %aln_seq;
open IN,"<$opts{fas}";
my $name="";
while (my $aline=<IN>){
	chomp $aline;
	if($aline=~/^>(\S+)/){
		$name=$1;
		next;
	}
	
	$aln_seq{$name} .=$aline;
}
close IN;

#reference
my $ref_acc=$opts{'ref'};
if(!defined $aln_seq{$ref_acc}){
	print "ERROR: Invalid reference sequence ID for -ref.\n";
	exit;
}

my @refseq=split //,$aln_seq{$ref_acc};
my $rawref=$aln_seq{$ref_acc};
$rawref=~s/-//g;
my $ref_length=length($rawref);

### Variant calling
my $gvf=$opts{"o"}."/mutations_result.gvf";
my $anno=$opts{"o"}."/mutations_stat.txt";
open GVF,">$gvf";
print GVF "$ref_acc\t$ref_acc\t\n";
open ANNO,">$anno";
print ANNO "#Ref_ID\tRef_length\tQuery_ID\tQuery_length\tSNP#\tInsertion#\tDeletion#\tIndel#\tSimilarity\n";
print ANNO "$ref_acc\t$ref_length\t$ref_acc\t$ref_length\t0\t0\t0\t0\t100\n";


foreach my $key (keys %aln_seq){
	my $Var_tag=0;
	my %sta_hash;
	my $P1=0;
	my $P2=0;
	my @pos;
	next if($key eq $ref_acc);
	my @queseq=split //,$aln_seq{$key};
	my @refseq_tmp=@refseq;
	my $MutationDetail="";
	my $query_seq = $aln_seq{$key};
	$query_seq=~s/-//g;
	my $query_length=length($query_seq);
	
	for (my $i=0;$i<=$#refseq_tmp ;$i++) {
		if($refseq_tmp[$i] eq "-" && (!$queseq[$i])){
			die "ERROR: The Sequence has some messy codes!\n";
			##new line: 2022-02-09 add.
		}else{
			if($refseq_tmp[$i] eq "-" && ($queseq[$i]=~/N/i|| $queseq[$i]=~/X/i)){
				my $NUMBER=$i;
				splice @refseq_tmp,$NUMBER,1;
				splice @queseq,$NUMBER,1;
				$i-=1;
				next;
			}
		}
	}
	for (my $i=0;$i<=$#refseq_tmp ;$i++) {
		if ($refseq_tmp[$i]=~/[ACGT]/i) {
			$P1++;
		}
		if ($queseq[$i]!~/\-/) {
			$P2++;
		}
		if ($queseq[$i]=~/N/ || $queseq[$i]=~/X/) {
			$queseq[$i]=$refseq_tmp[$i];
		}
		#next if($P1==0 || $P2==0);
		#last if(($P1==$reference_length && $P2<$query_length) || ($P1<$reference_length && $P2==$query_length));
		if ($refseq_tmp[$i] ne $queseq[$i]) {
			my $type="SNP";

			if ($refseq_tmp[$i] eq "-") {
				############Insertion############
				$type="Insertion";
				my ($refalle,$alle,$start);

				if ($P1==0) {
					$refalle=$refseq_tmp[0];
					$alle=$queseq[0];
				}else{
					$refalle=$refseq_tmp[$i-1];
					$alle=$queseq[$i-1];
				}
				if ($refalle ne $alle) {
					$type="Indel";
				}

				$start=$P1;
				my $n=0;
				for (my $j=$i;$j<=$#refseq_tmp ;$j++) {
					if($refseq_tmp[$j] ne "-"){
						if ($P1==0) {
							$refalle=$refseq_tmp[$j];
							$start=1;
							$alle.=$queseq[$j];
							$P1+=1;
							$i+=1;
						}
						last;
					}
					if ($P1==0) {
						if ($i!=$j) {
							$alle.=$queseq[$j];
						}
					}else{
						$alle.=$queseq[$j];
					}
					
					$n++;
				}
				$i+=$n-1;
				my $end=$P2-1;
				if ($P1==0) {
					$end=1;
				}
			
				push @pos,[$start,$start,$type];
				$MutationDetail.="$start($type:$refalle->$alle);";
				$Var_tag++;
				$P2=$P2+$n-1;
			}else{
				############Deltion############
				if ($queseq[$i] eq "-") {
					$type="Deletion";
					my ($refalle,$alle,$start);
					if($P2==0){
						$refalle=$refseq_tmp[0];
						$alle=$queseq[0];
						$start=1;
					}else{
						$refalle=$refseq_tmp[$i-1];
						$alle=$queseq[$i-1];
						if ($refalle ne $alle) {
							$type="Indel";
						}
						$start=$P1-1;
					}

					my $n=0;
					for (my $j=$i;$j<=$#queseq ;$j++) {
						last if($queseq[$j] ne "-");
						$refalle.=$refseq_tmp[$j];
						$n++;
					}
					$i+=$n-1;
					my $end=$P1+$n-1;
					
					push @pos,[$start,$end,$type];
					$MutationDetail.="$start($type:$refalle->$alle);";
					$Var_tag++;
					$P1=$end;
				}else{
					############SNP############
					my $k=$i+1;
					if ($k<=$#refseq_tmp) {
						next if(($refseq_tmp[$k] eq "-") || ($queseq[$k] eq "-"));
					}
					
					push @pos,[$P1,$P1,$type];
					$MutationDetail.="$P1($type:$refseq_tmp[$i]->$queseq[$i]);";
					$Var_tag++;
				}
			}
			$sta_hash{$type}++;
		}
	}

#print out
	my $TotalVarNo=0;
	####总突变统计
	if (!defined $sta_hash{"SNP"}) {
		$sta_hash{"SNP"}=0;
	}
	if (!defined $sta_hash{"Insertion"}) {
		$sta_hash{"Insertion"}=0;
	}
	if (!defined $sta_hash{"Deletion"}) {
		$sta_hash{"Deletion"}=0;
	}
	if (!defined $sta_hash{"Indel"}) {
		$sta_hash{"Indel"}=0;
	}

	foreach my $key (keys %sta_hash) {
		$TotalVarNo+=$sta_hash{$key};
	}
	my $r=1-2*$TotalVarNo/($ref_length+$query_length);
	$r=$r*100;
	$r=sprintf "%.2f",$r;

	####去掉头尾非SNP，统计中间突变
	if ($Var_tag>0) {
		@pos=sort {$a->[0]<=>$b->[0]} @pos;
		if (($pos[0][0]==1) && ($pos[0][2] ne "SNP")) {
			shift @pos;
		}
		my $pos_N=$#pos+1;
		if (($pos_N>=1) && ($pos[-1][1]==$ref_length) && ($pos[-1][2] ne "SNP")) {
			pop @pos;
		}
	}

	my $TV_middle=$#pos+1;

	my $ND=0;
	my $DL=0;
	my $Mid_var_length=0;
	my $Gap_num=0;
	for (my $i=0;$i<=$#pos ;$i++) {
		my $l=$pos[$i][1]-$pos[$i][0]+1;
		if ($pos[$i][2] eq "Deletion") {
			$ND++;
			if ($DL<$l) {
				$DL=$l;
			}
		}
		if ($pos[$i][2] ne "SNP") {
			$Gap_num++;
		}
		$Mid_var_length+=$l;
	}
	
	my $mutation_file="$key\t$key";
	if ($MutationDetail) {
		chop $MutationDetail;
		$mutation_file.="\t$MutationDetail";
	}
	$mutation_file.="\n";
	print GVF $mutation_file;
	
	print ANNO "$ref_acc\t$ref_length\t$key\t$query_length\t$sta_hash{\"SNP\"}\t$sta_hash{\"Insertion\"}\t$sta_hash{\"Deletion\"}\t$sta_hash{\"Indel\"}\t$r\n";
}


sub usage{
print <<"USAGE";
Version $version
Description: 
Usage:
$0 -fas -ref -o
options:
-fas input file, multiple aligned sequences
-ref input string, reference sequence id
-o output directory
-h help
USAGE
exit(1);
}


