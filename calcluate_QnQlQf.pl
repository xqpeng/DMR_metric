=cut
	function: 1. calculate the methyation difference and count the CpG number and the length in each DMR
				2.  Qn, Ql, and Qf under different threshold for a DMR set
	input: 1. a file containing formated DMR files of different methods on the same dataset(for each DMR file, the format is : chr	start_pos	end_pos	...)
			2. a merged methylation file of group A and group B
			3. 1 or 2 (used in calculating the methylation difference, e.g. 1 ---count based; 2 ---methylation level based;   ) 
	output: 1. a set of files corresponding to the formated DMR files, containing DMRs and the additional information, including the methyation difference and count the CpG number and the length of each DMR
			2. a file containing Qn, Ql, and Qf of different methods under different threshold
=cut

use strict;
use locale; 
use POSIX qw(log10);
my $DMRfile_list=shift;
my $merged_meth_groups=shift;
my $meth_dif_method=shift; #choose meth_difference calculation method.

my %Ql=();
my %Qn=();
my %Qf=();
my %DMR_th=();
my %probe_methA=();
my %probe_ummeA=();
my %probe_levelA=();
my %probe_sampleA=();
my %probe_methB=();
my %probe_ummeB=();
my %probe_levelB=();
my %probe_sampleB=();
open(RF,$merged_meth_groups)or die $!;
while(my $line=<RF>){
	chomp $line;
	my @temp=split /\t/,$line;
	$probe_sampleA{$temp[0]}{$temp[1]}=$temp[2];
	$probe_levelA{$temp[0]}{$temp[1]}=$temp[3];
	$probe_methA{$temp[0]}{$temp[1]}=$temp[4];
	$probe_ummeA{$temp[0]}{$temp[1]}=$temp[5];
	$probe_sampleB{$temp[0]}{$temp[1]}=$temp[6];
	$probe_levelB{$temp[0]}{$temp[1]}=$temp[7];
	$probe_methB{$temp[0]}{$temp[1]}=$temp[8];
	$probe_ummeB{$temp[0]}{$temp[1]}=$temp[9];
}
close(RF);

open(LS,$DMRfile_list)or die $!;
my %DMRset=();
while(my $DMRfile=<LS>){
	chomp $DMRfile;
	my $output="Meth_Dif_".$DMRfile;
	my $chr="";
	print "process...$DMRfile\n";
	%DMRset=();
	open(IN, $DMRfile)or die $!;
	open(OUT,">$output")or die $!;
	while(my $line=<IN>){
		chomp $line;
		my @temp=split /\t/,$line;
		if($chr ne $temp[0]){
			if($chr ne ""){
				my $output_string=&process_chr_DMR($chr);
				print OUT $output_string;
				%DMRset=();
				
			}
			$chr=$temp[0];
		}
		$DMRset{$temp[1]}=$temp[2];
	}
	my $output_string=&process_chr_DMR($chr);
	print OUT $output_string;
	close(IN);
	close(OUT);

	
}
close(LS);
open(LS,$DMRfile_list)or die $!;
while(my $DMRfile=<LS>){
	chomp $DMRfile;
	&interval_statistical($DMRfile,$meth_dif_method);
}
close(LS);
my @method_results=keys(%Qn);
my %method_name=();
for my $each(@method_results){
	my @temp=split /_/,$each;
		$method_name{$temp[1]}=$each;
	}
my @method=sort keys(%method_name);
my $title_line=join("\t",@method);
&print_result($meth_dif_method);



sub print_result{
	my $method_i=shift;
	print "original: $title_line\n";
	my $Qn_results="Qn\t".$title_line."\n";
	my $Ql_results="Ql\t".$title_line."\n";
	my $Qf_results="Qf\t".$title_line."\n";
	my $DMR_results="Ratio\t".$title_line."\n";
	open(OUT1,">>QnQlQf_original.txt")or die $!;
	for(my $i=0;$i<1;$i+=0.1){
		$Qn_results.=$i;
		$Ql_results.=$i;
		$Qf_results.=$i;
		$DMR_results.=$i;
		for my $each_method(@method){
			my $each=$method_name{$each_method};			
			$Qn_results.="\t".$Qn{$each}{$i};
			$Ql_results.="\t".$Ql{$each}{$i};
			$Qf_results.="\t".$Qf{$each}{$i};
			$DMR_results.="\t".$DMR_th{$each}{$i};
		}
		$Qn_results.="\n";
		$Ql_results.="\n";
		$Qf_results.="\n";
		$DMR_results.="\n";
	}
	print $Qn_results;
	print OUT1 "\n\n\n".$DMRfile_list.":Qn, Ql, Qf\n";
	if($method_i==1){
		print OUT1 "1 ---count based:\n"
	}else{
		print OUT1 "2 ---methylation level based:\n";
	}
	print OUT1 $merged_meth_groups."\n";
	print OUT1 $Qn_results;
	print OUT1 $Ql_results;
	print OUT1 $Qf_results;
	print OUT1 $DMR_results;
	close(OUT1);
}

sub process_chr_DMR{
	my $curr_chr=shift;
	my @probes=sort{$a<=>$b}keys(%{$probe_sampleA{$curr_chr}});
	
	my @DMR_start=sort{$a<=>$b}keys(%DMRset);
	#print "$curr_chr\t".@DMR_start."\t".@probes."\n";
	my $cur=0;
	my $return_inf="";
	my %DMR_cpg=();
	my %DMR_length=();
	my %DMR_sample_A=();
	my %DMR_sample_B=();
	my %DMR_difmeth_1=();
	my %DMR_difmeth_2=();
	my %DMR_difmeth_1_meth_A=();  #count based differetial methylation
	my %DMR_difmeth_1_unmeth_A=(); 
	my %DMR_difmeth_1_meth_B=();  #count based differetial methylation
	my %DMR_difmeth_1_unmeth_B=(); 
	my %DMR_difmeth_2_A=();  #methylation level based differential methylation
	my %DMR_difmeth_2_B=();  #methylation level based differential methylation
	for my $start(@DMR_start){		
		$DMR_cpg{$start}=0;
		$DMR_length{$start}=0;
		$DMR_sample_A{$start}=0;
		$DMR_sample_B{$start}=0;
		$DMR_difmeth_1_meth_A{$start}=0;
		$DMR_difmeth_1_unmeth_A{$start}=0;
		$DMR_difmeth_2_A{$start}=0;
		$DMR_difmeth_1_meth_B{$start}=0;
		$DMR_difmeth_1_unmeth_B{$start}=0;
		$DMR_difmeth_2_B{$start}=0;
		$DMR_difmeth_1{$start}=0;
		$DMR_difmeth_2{$start}=0;
		my %cpg_interval=();
		my @dif_interval=();
		for (my $i=$cur;$i<@probes;$i++){
			if($probes[$i]>$DMRset{$start}){
				$cur=$i;
				#print "$curr_chr\t$probes[$i]\t$DMRset{$start}\n";
				last;
			}elsif($probes[$i]<$start){
				next;
			}else{
				
				$DMR_cpg{$start}++;
				$DMR_sample_A{$start}=$DMR_sample_A{$start}+$probe_sampleA{$curr_chr}{$probes[$i]};
				$DMR_sample_B{$start}=$DMR_sample_B{$start}+$probe_sampleB{$curr_chr}{$probes[$i]};
				if($probe_sampleA{$curr_chr}{$probes[$i]}!=0){
					$DMR_difmeth_1_meth_A{$start}=$DMR_difmeth_1_meth_A{$start}+$probe_methA{$curr_chr}{$probes[$i]};
					$DMR_difmeth_1_unmeth_A{$start}=$DMR_difmeth_1_unmeth_A{$start}+$probe_ummeA{$curr_chr}{$probes[$i]};
					$DMR_difmeth_2_A{$start}=$DMR_difmeth_2_A{$start}+$probe_sampleA{$curr_chr}{$probes[$i]}*$probe_levelA{$curr_chr}{$probes[$i]};
				}
				if($probe_sampleB{$curr_chr}{$probes[$i]}!=0){
					$DMR_difmeth_1_meth_B{$start}=$DMR_difmeth_1_meth_B{$start}+$probe_methB{$curr_chr}{$probes[$i]};
					$DMR_difmeth_1_unmeth_B{$start}=$DMR_difmeth_1_unmeth_B{$start}+$probe_ummeB{$curr_chr}{$probes[$i]};
					$DMR_difmeth_2_B{$start}=$DMR_difmeth_2_B{$start}+$probe_sampleB{$curr_chr}{$probes[$i]}*$probe_levelB{$curr_chr}{$probes[$i]};
				}
				
				if($probe_sampleA{$curr_chr}{$probes[$i]}!=0 && $probe_sampleB{$curr_chr}{$probes[$i]}!=0){
					my $pos_dif=abs($probe_levelA{$curr_chr}{$probes[$i]}-$probe_levelB{$curr_chr}{$probes[$i]});
					push @dif_interval,$pos_dif;

					$pos_dif = int($pos_dif * 10);
					if(exists($cpg_interval{$pos_dif/10})){
						$cpg_interval{$pos_dif/10}+= 1;
					}else{
						$cpg_interval{$pos_dif/10}= 1;
					}
				}
			}
		}
		if(exists($cpg_interval{1})){
			if(exists($cpg_interval{0.9})){
				$cpg_interval{0.9}+=$cpg_interval{1};
			}else{
				$cpg_interval{0.9}=$cpg_interval{1};
			}
		}
		my @interval=keys(%cpg_interval);
		
		my $total_count = 0;
		my $sum_ratio = 0;
		for(my $i = 0; $i <1; $i = $i + 0.1){
			if(!exists($cpg_interval{$i})){
				$cpg_interval{$i} = 0;
			}
			$sum_ratio=$sum_ratio+($i+0.05);
			$total_count=$total_count+($i+0.05)*$cpg_interval{$i};
		}
		
		my $DMR_rank=$total_count/$sum_ratio;
		
		if($DMR_sample_A{$start}>0 && $DMR_sample_B{$start}>0){
			#print $curr_chr."\t".$start."\t".$DMRset{$start}."\t".($DMRset{$start}-$start+1)."\t".$DMR_cpg{$start}."\t".$DMR_sample_A{$start}."\t".$DMR_sample_B{$start}."\t".($DMR_difmeth_1_meth_A{$start}+$DMR_difmeth_1_unmeth_A{$start})."\t".($DMR_difmeth_1_meth_B{$start}+$DMR_difmeth_1_unmeth_B{$start})."\n";
			$DMR_difmeth_1{$start}=abs($DMR_difmeth_1_meth_A{$start}/($DMR_difmeth_1_meth_A{$start}+$DMR_difmeth_1_unmeth_A{$start})-$DMR_difmeth_1_meth_B{$start}/($DMR_difmeth_1_meth_B{$start}+$DMR_difmeth_1_unmeth_B{$start}));
			$DMR_difmeth_2{$start}=abs($DMR_difmeth_2_A{$start}/$DMR_sample_A{$start}-$DMR_difmeth_2_B{$start}/$DMR_sample_B{$start});
		}else{
		#	print $curr_chr."\t".$start."\t".$DMRset{$start}."\t".($DMRset{$start}-$start+1)."\t".$DMR_cpg{$start}."\t".$DMR_sample_A{$start}."\t".$DMR_sample_B{$start}."\t".($DMR_difmeth_1_meth_A{$start}+$DMR_difmeth_1_unmeth_A{$start})."\t".($DMR_difmeth_1_meth_B{$start}+$DMR_difmeth_1_unmeth_B{$start})."\n";
			$DMR_difmeth_1{$start}=0;
			$DMR_difmeth_2{$start}=0;
		}
		$return_inf.=$curr_chr."\t".$start."\t".$DMRset{$start}."\t".($DMRset{$start}-$start+1)."\t".$DMR_cpg{$start}."\t".$DMR_rank."\t".$DMR_difmeth_2{$start}."\t".$DMR_difmeth_1{$start}."\n";
	}
	return $return_inf;
}


sub interval_statistical{
	my $file=shift;
	my $method=shift;#choose meth_difference calculation method.
	
	my $input_file="Meth_Dif_".$file;
	my $minum_cpg=5; #default=5
	my $maxium_length=10000; #default=10000
	my %count_interval=();
	my %length_interval=();
	
	open(IN,$input_file)or die $!;
	while(my $line= <IN>){
		chomp $line;
		my @temp = split /\t/,$line;
		if($temp[4] <$minum_cpg || $temp[3] >=$maxium_length){
			next;
		}
		
		my $c = int($temp[-$method] * 10);
		#if($temp[-$method]>0.5){
		#	print "$temp[-$method]\t$c\t".($c/10)."\n";
		#}
		if(exists($count_interval{$c/10})){
			$length_interval{$c/10}+=$temp[3];
			$count_interval{$c/10}+=$temp[4];
		}else{
			$length_interval{$c/10}=$temp[3];
			$count_interval{$c/10}=$temp[4];
		}
		if(!exists($DMR_th{$file}{$c/10})){
			$DMR_th{$file}{$c/10}=1;
		}else{
			$DMR_th{$file}{$c/10}++;
		}
	}
	close(IN);
	if(exists($count_interval{1})){
		if(exists($count_interval{0.9})){
			$count_interval{0.9}+=$count_interval{1};
			$length_interval{0.9}+=$length_interval{1};
			$DMR_th{$file}{0.9}+=$DMR_th{$file}{1};
		}else{
			$count_interval{0.9}=$count_interval{1};
			$length_interval{0.9}=$length_interval{1};
			$DMR_th{$file}{0.9}=$DMR_th{$file}{1};
		}
	}
	my $total_count=0;
	my $total_length=0;
	my $sum_ratio=0;
	my @i_value=sort{$a<=>$b}keys(%count_interval);
	print "@i_value\n";
	for(my $i=0.9;$i>=0;$i=$i-0.1){
		if(!exists($count_interval{$i})){
			$count_interval{$i}=0;
			$length_interval{$i}=0;
			$DMR_th{$file}{$i}=0;
		}
		$sum_ratio=$sum_ratio+($i+0.05);
		$total_count=$total_count+($i+0.05)*$count_interval{$i};
		$total_length=$total_length+($i+0.05)*$length_interval{$i};
		if($total_count>0){
			$Qn{$file}{$i}=$total_count/$sum_ratio;
			$Ql{$file}{$i}=$total_length/$sum_ratio;
		}else{
			$Qn{$file}{$i}=0;
			$Ql{$file}{$i}=0;
		}
		if(($Qn{$file}{$i}+$Ql{$file}{$i})==0){
			$Qf{$file}{$i}=0;
		}else{
			$Qf{$file}{$i}=2*$Ql{$file}{$i}*$Qn{$file}{$i}/($Qn{$file}{$i}+$Ql{$file}{$i});
		}
		$DMR_th{$file}{$i}=$DMR_th{$file}{$i}+$DMR_th{$file}{$i+0.1};
		if($DMR_th{$file}{$i}>0){
	#		print "$file\t$i\t$DMR_th{$file}{$i}\t$total_count\t$total_length\t$sum_ratio\t$Qn{$file}{$i}\t$Ql{$file}{$i}\n";
		}
		
		
	}
}

