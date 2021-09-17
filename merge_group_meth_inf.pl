=cut
	function: Merge the methylation report files of samples in each group into one methylation report file;
			It provides two ways of merging for each cytosine: 1. count based: merging the numbers of methylated and unmethylated cytosines from each sample;
											2. methylation level based: calculating the average methylation level based on the methylation levels from samples.
	Input: 1. A file containing the methylation report files of samples with group id
			For example: methylation_bedGraph_sample1 A
						methylation_bedGraph_sample2 A
						methylation_bedGraph_sample3 B
						methylation_bedGraph_sample4 B
						...
						(methylation_bedGraph_sample1 is a bedGraph/coverage output file, 
						the format is :chr	pos	strand	meth	unmeth)
			2. the column id to specify the meth counts
			3. the column id to specify the unmeth counts
	Output: merged_methylation_bedGraph_A_B
			
=cut


use strict;



my $list=shift;
my $meth_col=shift;
my $unmeth_col=shift;
#my $meth_col=3;
#my $unmeth_col=4;
my $filter_depth=1;
my %samp=();
my %all_chr=();
my $out="merged_v9_".$list;
open(OUTT,">$out")or die $!;		
print OUTT "Chr\tpos\tsample_numA\tave_meth_level_A\tmethA\tunmethA\tsamplenum_B\tave_meth_level_B\tmethB\tunmethB\n";

open(LI,$list)or die $!;

while(my $file=<LI>){
	chomp $file;	
	my @samp_inf=split /\t/,$file;
	$samp{$samp_inf[1]}{$samp_inf[0]}=1;
	
}
close(LI);

my @class=sort{$a<=>$b}keys(%samp);
print "Groups @class\n";

my %probe_meth=();
my %probe_umme=();
my %probe_level=();
my %probe_sample=();
my %probe_infA=();
my %probe_infB=();
my %probe=();
my $N=0;
for my $type1(@class){
	
	$N++;
	my @sam_file=keys(%{$samp{$type1}});
	%probe_meth=();
	%probe_umme=();
	%probe_level=();
	%probe_sample=();
	for my $file(@sam_file){
		print "processing....$file\n";
		
		open(IN,$file)or  die;
		 
		while(my $line=<IN>){
			chomp $line;
			my @temp=split /\t/,$line;
			if($temp[$meth_col]+$temp[$unmeth_col]>=$filter_depth){
				$probe{$temp[0]}{$temp[1]}=1;					
				if(exists($probe_sample{$temp[0]}{$temp[1]})){
					$probe_sample{$temp[0]}{$temp[1]}++;
					$probe_meth{$temp[0]}{$temp[1]}+=$temp[$meth_col];
					$probe_umme{$temp[0]}{$temp[1]}+=$temp[$unmeth_col];
					$probe_level{$temp[0]}{$temp[1]}+=$temp[$meth_col]/($temp[$meth_col]+$temp[$unmeth_col]);
				}else{
					$probe_sample{$temp[0]}{$temp[1]}=1;
					$probe_meth{$temp[0]}{$temp[1]}=$temp[$meth_col];
					$probe_umme{$temp[0]}{$temp[1]}=$temp[$unmeth_col];
					$probe_level{$temp[0]}{$temp[1]}=$temp[$meth_col]/($temp[$meth_col]+$temp[$unmeth_col])
				}
			}
			
		}
		close(IN);
	}
	my @chromo=sort{$a<=>$b}keys(%probe_sample);
	for my $chr(@chromo){
		my @probe_key=sort{$a<=>$b}keys(%{$probe_sample{$chr}});
		
		if($N==1){
			for my $prb(@probe_key){
				$probe_infA{$chr}{$prb}=$probe_sample{$chr}{$prb}."\t".($probe_level{$chr}{$prb}/$probe_sample{$chr}{$prb})."\t".$probe_meth{$chr}{$prb}."\t".$probe_umme{$chr}{$prb};
					
			}
		}else{
			for my $prb(@probe_key){
				$probe_infB{$chr}{$prb}=$probe_sample{$chr}{$prb}."\t".($probe_level{$chr}{$prb}/$probe_sample{$chr}{$prb})."\t".$probe_meth{$chr}{$prb}."\t".$probe_umme{$chr}{$prb};
			}
		}
	
	}
	
}
my @chromo=sort{$a<=>$b}keys(%probe);

for my $chr(@chromo){
	my @probe_key=sort{$a<=>$b}keys(%{$probe{$chr}});
	for my $prb(@probe_key){
		if(exists($probe_infA{$chr}{$prb}) && exists($probe_infB{$chr}{$prb})){
			print OUTT $chr."\t".$prb."\t".$probe_infA{$chr}{$prb}."\t".$probe_infB{$chr}{$prb}."\n";
		}elsif(exists($probe_infA{$chr}{$prb})){
			print OUTT $chr."\t".$prb."\t".$probe_infA{$chr}{$prb}."\t0\tNA\tNA\tNA\n";
		}else{
			print OUTT $chr."\t".$prb."\t0\tNA\tNA\tNA\t".$probe_infB{$chr}{$prb}."\n";
		}
	}
}

close(OUTT);


