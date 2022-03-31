#!/usr/bin/env perl

use Modern::Perl "2011";
use autodie;
use Smart::Comments;
use Tie::IxHash;
use Storable;
use Data::Dumper;

unless (@ARGV == 5) {
    die <<"EOT";
Usage: $0 <meme.txt> <number> <date> <suffix> <indir>
This script creates the matrices associated with MEME.txt files.
It requires 3 arguments:
(1) meme.txt file (full path)
(2) which MEME motif was selected (1, 2, 3, or 4)
(3) the date
(4) the suffix
(5) indir (ex. "~/Documents/bioinfo/COG_directory/OG-")
It creates a matrix Fasta file for the given transcription factor COG.
Used in a loop it can create several matrices.
Example: 
for path in `cat meme_links_180724.txt`
 do meme_file=`echo \$path | cut -d" " -f1`
 number=`echo \$path | cut -d" " -f2`
 perl /Users/sinaeda/Documents/bioinfo/perl_scripts/extract_matrix.pl \$meme_file \$number 180724 \$indir; done
done

EOT
}

# Used in for loop

my $meme_file = shift; 	# ~/Documents/bioinfo/ProteinOrtho/COG90/LacI/SCO0456_UPS/OG-SCO0456_meme/CDS300_zoops_medium/meme.txt
my $motif_num = shift; # LacI_COG.dat 	
my $date = shift;
my $suffix = shift;
my $indir = shift;
my @family = $meme_file =~ /.*\/COG90\/(.*?)\//; # don't use this line pfiew
my @temp = split "/", $meme_file;
my @species = $meme_file =~ /.*\/OG-(.*)_UPS\/.*/;
my $laci_file = $indir . $species[0] . ".fasta";
my $outfile = $species[0]."_motif".$motif_num."_".$suffix."_$date.fasta";

my $num_rgx = "Motif .+ MEME-$motif_num in BLOCKS format";

# Open different in and outfiles
open my $in, '<', $meme_file;
open my $laci, '<', $laci_file;
open my $out, '>', $outfile;

# Create array where @fields[1] is the sequence number (e.g. 1_1111)
LINE:
#my $lines;
my @tf_ids;
tie my %motif_for, 'Tie::IxHash';
tie my %strand_for, 'Tie::IxHash';
my @strands;
while (<$in>) {
	chomp; 
	# Read what is between the num_rgx, and the end of the motifs (//), line by line
	if (/Motif .* MEME-$motif_num sites sorted by position p-value/../Motif .* MEME-$motif_num block diagrams/) {
		next if /Motif .* MEME-$motif_num sites sorted by position p-value/ or /Motif .* MEME-$motif_num block diagrams/ or /^-.*/ or /^Sequence.*/;
		my $line = $_;
		my @fields = split ' ', $line;
		my @id = $fields[0] =~ m/(\d+\_\d+)\(\-.*/xms if defined $fields[0];
		my $id = ">" . $id[0] if defined $id[0];
		my $strand = $fields[1] if defined $fields[1];
		
		push @{ $strand_for{$id} }, \$strand if defined $id;
	}
	if (/$num_rgx/../^\/\//) {
		# Don't read line if it is uninteresting (title line, line with -----, etc)
		next if /$num_rgx/ or /^\/\// or /^-.*/ or /^BL.*/; 
		my $line = $_ ;
		my @fields = split(/\(|\)/, $line) ;
		# tf_id is going to be the key in the hash %motif_for
		my $tf_id = ">" . $fields[0] ;
		
		# change $fields[0] for the pattern searching in map, to exclude any ambiguity 
		# (ex: 1_11 might match 1_110, whereas "1_11 lcl" will only match the same start)
		$fields[0] = ">" . $fields[0] . " lcl" ;
		push @tf_ids, $fields[0];

		# create hash with tf_id (key) and the rest of the info (motif, location, etc) as values inside array (new @fields)
		@fields = @fields[1,3,4];
		push @{ $motif_for{$tf_id} }, \@fields if defined $fields[0];
		## $tf_id
		## %motif_for
	}
}
## %strand_for
## %motif_for


# Create array of LacI_COG.dat
my @laci_tfs;
while (<$laci>) {
	chomp;
	if (/^>/) { push @laci_tfs, $_ }
}
## @laci_tfs

# Find the description of our LacIs
my @laci_descr = map {	my $tf_id = $_ ;
						grep { $_ =~ $tf_id } @laci_tfs;
				} @tf_ids;			
## @laci_descr

# Add to these descriptions a UPS and search for the corresponding motif and info
foreach(@laci_descr) {
	my @split_descr = split(" ", $_);
	my @locus_tag = grep (/\[locus_tag/, @split_descr);
	# $index should be matching keys of %motif_for
	my $index = $split_descr[0];
	my @prot_id = grep (/protein_id/, @split_descr);
	my @loc = grep (/location/, @split_descr);
	my @location = $loc[0] =~ /(\d+)\.\.(\d+)/;
	## $index
	## @loc
	## @location
	
	# Because we remove a key/value pair each time a key is used, we have to set up a condition checking 
	# that hash{$index} is defined or has already been used and removed from the hash
	if (defined @{$motif_for{$index}}[0]) {
		my @motif_info = @{$motif_for{$index}} ;
		my $n_motifs = scalar( @motif_info ) - 1 ;
		## $index
		## @motif_info
		## $n_motifs
		
		# OK for UPS_length because for the two motifs the UPS length is the same
		my $UPS_info = @{ $motif_info[0] }[0];
		my @UPS_length = split ',', $UPS_info;
		#my $UPS_length = $UPS_length[2];
		## $UPS_length
		## $UPS_info
		
		my @motif_start;
		my @motifs;
		for my $i (0..$n_motifs) {
			my @start = split ' ', @{ $motif_info[$i] }[1];
			push @motif_start, $start[0];
			my @motif = split ' ', @{ $motif_info[$i] }[2];
			push @motifs, $motif[0];
			
			# Strand 
			my $strand = @{ $strand_for{$index}}[$i];
			# $strand is a reference to the strand, so in the say expression, I dereference it by adding a $ --> $$strand
			
			if ($loc[0] =~ /.*complement.*/) {
				
				my $coordinates = ($location[1] + 1) . ".." . ($location[1] + $UPS_length[2]);
				## $coordinates
				say {$out} $index . " UPS " . $locus_tag[0] . " " . $prot_id[0] . " [length = ($UPS_info)] [motif_start = $motif_start[$i]] [coordinates = $coordinates] [".$$strand." strand]\n" . $motifs[$i];
			}
			else {
				my $coordinates = ($location[0] - $UPS_length[2]) . ".." . ($location[0] - 1);
				## $coordinates
				say {$out} $index . " UPS " . $locus_tag[0] . " " . $prot_id[0] . " [length = ($UPS_info)] [motif_start = $motif_start[$i]] [coordinates = $coordinates] [".$$strand." strand]\n" . $motifs[$i];
			}
		}
		## @motif_start
		## @motifs
	}
	# When one index has been used, remove the correspond key/value pair from hash so there will be no duplicates in the outfile
	delete $motif_for{$index};
}










