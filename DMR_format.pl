=cut
	Function: format and sort the DMR set predicted by different DMR detection tools
	Input: a file containing the names of DMR files of different methods on the same dataset
	Output: formated DMR files 
=cut 

use strict;
my $DMRfile_list=shift;
open(LS,$DMRfile_list)or die $!;
my %DMRset=();
while(my $line=<LS>){
	chomp $line;
	my $DMRfile;
	my $title_flag;
	my $chr_col;
	my $start_col;
	my $end_col;
	
	if($line=~/\s/){
		my @filed=split /\s/,$line;
		$DMRfile=$filed[0];
		$title_flag=$filed[1];
		$chr_col=$filed[2];
		$start_col=$filed[3];
		$end_col=$start_col+1;
	}else{
		$DMRfile=$line;
		$title_flag=1;
		$chr_col=0;
		$start_col=1;
		$end_col=2;
	}
	
	my $output="sort_".$DMRfile;
	%DMRset=();
	open(IN, $DMRfile)or die $!;
	open(OUT,">$output")or die $!;
	if($title_flag==1){
	#	my $title=<IN>;
	}	
	while(my $line=<IN>){
		chomp $line;
		my @temp=split /\t/,$line;
		if($temp[$chr_col]=~/^chr[\d+|X|Y|M]/){
			$DMRset{$temp[$chr_col]}{$temp[$start_col]}=$temp[$chr_col]."\t".$temp[$start_col]."\t".$temp[$end_col];
		}elsif($temp[$chr_col]=~/^[\d+|X|Y|M]/){
			$DMRset{$temp[$chr_col]}{$temp[$start_col]}="chr".$temp[$chr_col]."\t".$temp[$start_col]."\t".$temp[$end_col];
		}else{
		#	print "The title of ".$DMRfile.": $line\n";
			next;
		}
		
		for(my $i=0;$i<=@temp;$i++){
			if($i!=$chr_col && $i!=$start_col && $i!=$end_col){
			$DMRset{$temp[$chr_col]}{$temp[$start_col]}.="\t".$temp[$i];
			}
		}
		
	}
	my @chr=sort{$a<=>$b}keys(%DMRset);
	my $output_string="";
	for my $each_chr(@chr){
		my @DMR_start=sort{$a<=>$b}keys(%{$DMRset{$each_chr}});
		for my $pos(@DMR_start){
			$output_string.=$DMRset{$each_chr}{$pos}."\n";
		}
	}
	print OUT $output_string;
	close(IN);
	close(OUT);
	
}
close(LS);
