#!/usr/bin/env perl
use strict;

#### Extract SNPs from multiple alignment ####
### Return  position and sample in which SNP is unique ####

# Read input and output file names from command line
#my ($ref_file, $output_file5) = @ARGV;
if ($#ARGV != 4)
{
    print "Input and output files should be specified, as well as the name of the sequences being aligned.\n";
}
my ($input_file, $output_file, $seed1, $seed2, $seed3) = @ARGV;

## Change file names & path
#my $input_file  = "/Users/erika/Data/sftp_VSC_stg00113/NOVOLoci/FromNicolas/contigs_61373_61335_renamed.aln";
#my $output_file  = "/Users/erika/Data/sftp_VSC_stg00113/NOVOLoci/FromNicolas/contigs_61373_61335_SNPs_2.txt";

open(INPUT, $input_file) or die "Can't open file $input_file, $!\n";
open(OUTPUT, ">" .$output_file) or die "Can't open file $output_file, $!\n";

my %seed;

my $pos_ref = "patient";
## Name of index + parents
$seed{$seed1} = undef;
$seed{$seed2} = undef;
$seed{$seed3} = undef;
#$seed{"61335_LCRC_HAP161335_LCRC_HAP1"} = undef;
#$seed{"61335_LCRC_HAP2"} = undef;
#$seed{"61373_LCRD_HAP1"} = undef;


#merge mafft lines-------------------------------------------------------------------------------             
my $g = '0';
my $query_line = "";
my %subject_list;
undef %subject_list;
my $consensus_total = "";
my %gaps_id;
undef %gaps_id;
my %length_id;
undef %length_id;

INPUT: while (my $line2 = <INPUT>)
{                                                     
    chomp($line2);
    if ($g > 2)
    {
        my @blast_result_tmp = split /\s+/, $line2;
        # Store sequenced of all alignments in hash
        if (exists($seed{$blast_result_tmp[0]}) && $blast_result_tmp[0] ne "")
        {
            my $subject_tmp = $subject_list{$blast_result_tmp[0]};
            $subject_list{$blast_result_tmp[0]} = $subject_tmp.$blast_result_tmp[1];
            $query_line = "yes";
        }
	# Store consensus in hash
        elsif ($query_line eq "yes")
        {
            my $consensus = substr $line2, 16, 60;
            $consensus_total .= $consensus;
        }
    }
    $g++
}
my $v = '0';
my $prev_pos = "";
PRINT: while ($v < length($consensus_total)-1)
{
   my $pos = substr $consensus_total, $v, 1;
   my $pos_next = substr $consensus_total, $v+1, 1;
   my $print_text = "";
   my %nucs;
   undef %nucs;

   my %list_subject;
   my %count_base;
   my %count_non_empty;
   my %count_duplicate;
   my $position;
   
   # Check positions with variants (some positions with variants are not indicated by ".")
   if ($pos ne "*" && $prev_pos eq "*" && $pos_next eq "*")
   #if ($pos eq "." && $prev_pos eq "*" && $pos_next eq "*")
   {
       undef %list_subject;
       undef %count_base;
       #undef %count_non_empty;
       
       print $pos."\n";
       my $seq_tmp = substr $subject_list{$pos_ref}, 0, $v;
       my $gaps = $seq_tmp =~ tr/-/-/;
       $position = $v+1-$gaps;
       $count_non_empty{$position} = 0;
       $count_duplicate{$position} = 0;
       
       # Check sequence of all samples in alignment
       foreach my $nuc_subject (keys %subject_list)
       {				
	   my $pos_sub = substr $subject_list{$nuc_subject}, $v, 1;
	   print "CHECK IN\t".$position."\t".$nuc_subject."\t".$count_base{$pos_sub}."\t".$list_subject{$pos_sub}."\n";
	   
	   if ($pos_sub ne "-" && $pos_sub ne "n" && $pos_sub ne "N" && $pos_sub ne "")
	   {
	       $count_non_empty{$position}++;
	       print $nuc_subject."\t".$pos_sub."\n";
	       if (exists $count_base{$pos_sub})
	       {
		   $count_base{$pos_sub}++;
		   $list_subject{$pos_sub} .= $nuc_subject;
		   $count_duplicate{$position}++;
	       }
	       else
	       {
		   $count_base{$pos_sub} = 1;
		   $list_subject{$pos_sub} = $nuc_subject;
	       }
	       print "CHECK OUT\t".$pos_sub."\t".$nuc_subject."\t".$count_base{$pos_sub}."\t".$list_subject{$pos_sub}."\n";
	   }   
	   # Store in print_text base with variant, not in first sequence
	   #if (exists($nucs{$pos_sub}) && $pos_sub ne "-" && $pos_sub ne "n" && $pos_sub ne "N")
	   #{
	   #    $print_text = $pos_sub;
	   #}
	   # Only store in nucs hash name of 1st sample with variant
	   #elsif ($pos_sub ne "-" && $pos_sub ne "n" && $pos_sub ne "N")
	   #{
	   #    $nucs{$pos_sub} = $nuc_subject;
	   #}
	   #else
	   #{
	   #    $v++;
	   #    next PRINT;
	   #}
	   }

       if ($count_non_empty{$position} > 1 && $count_duplicate{$position} >= 1)
       { 
	   foreach my $key (keys %count_base)
	   {
	       #if ($key 
	       if ($count_base{$key} == 1)
	       {
		   $nucs{$key} = $list_subject{$key};
		   $print_text = $key;
		   print OUTPUT $position."\t";
		   print OUTPUT $list_subject{$key}."\n";
	       }
	   }
       } 
       print "FINAL\t".$print_text."\t".$nucs{$print_text}."\n";
       #if ($print_text ne "" && $print_text ne "n" && $print_text ne "N" && $print_text ne "")
       #{
	#   if ($count_non_empty{$position} > 1)
	#   { 
	#   foreach my $nuc_check (keys %nucs)
	#   {
	       #if ($nuc_check ne $print_text)
	       #{
		   #my $gg = $v+1;
		   #my $seq_tmp = substr $subject_list{$pos_ref}, 0, $v;
		   #my $gaps = $seq_tmp =~ tr/-/-/;
		   #my $position = $v+1-$gaps;
	#	   print OUTPUT $position."\t";
	#	   print OUTPUT $nucs{$nuc_check}."\n";
		   #}
	 #  }	   
	 #  }
       #}
   }
   $prev_pos = $pos;
   $v++;
}

close INPUT;
close OUTPUT;
