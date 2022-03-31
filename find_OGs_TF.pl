#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;
use Smart::Comments;
use Tie::IxHash;
use List::AllUtils qw(uniq);
use Storable;

my $begin = time;

# Create tf_ids, which opens and then closes the dat file for the TF
# Only once
my $indir = shift; # directory with the COG.fasta files
my $dat = shift;
my $outfile = shift;

# Open directory and get all the files in it (perl has no 'argument list too long' error)
opendir(DIR, $indir);
my @infiles = grep(/OG-.*\.fasta/, readdir(DIR));
closedir(DIR);

open my $dh, '<', $dat;
my @tf_ids;
while (my $l = <$dh>) {
	chomp $l;
	push @tf_ids, $l;
}
close $dh;


# Create hash for COGs and their corresponding gene IDs
tie my %ids_for, 'Tie::IxHash';

if (! -e join($indir, 'COG_hash.ref')) {
	for my $infile (@infiles) {
		## $infile
		%ids_for = (%ids_for, read_COG_fasta($infile));
	}
	## %ids_for
	store \%ids_for, join($indir, 'COG_hash.ref');
}

my $id_ref = retrieve(join($indir, 'COG_hash.ref'));
%ids_for = %$id_ref;
## %ids_for


# Go through $dat to look for the same ">number_id lcl" and keep the COGs
# We have to get something like this:
# ./OG-XNR_RS03065.fasta:>50_5124 lcl [protein=transcriptional regulator]

my @keys = keys %ids_for;
## @keys

open my $out, '>', $outfile; 

my @goods = map { 	my @vals = @{$ids_for{$_}};
					my @matches = map { my $tf_id = $_;
										grep {$_ =~ $tf_id} @vals } @tf_ids;
					my $num_match = scalar(@matches);
					## $num_match
					for my $match ( @matches[0..($num_match-1)] ){
						say {$out} "./" . $_ . ".fasta:" . $match if (scalar(@matches) > 0);
					} 
				} @keys;

## @goods

sub read_COG_fasta {			# recycled "as is" from xxl_xlate.pl
	my $infile = shift;
	my ($cog, $suffix) = $infile =~ m/(OG-.*)\.fasta/xms;
	open my $in, '<', $infile;
	my $seq_id;
	my @ids;
	tie my %ids_for, 'Tie::IxHash';
	LINE:
	while (my $line = <$in>) {
		chomp $line; 
		if ($line =~ m/^>/xms) {
			my @id = $line =~ m/(>\d+_\d+\s+lcl).*protein=(.*)]\s+\[protein_id/xms;
			## @id
			$seq_id = $id[0]." [protein=".$id[1]."]";
			## $seq_id
		}
		push @ids, $seq_id;
	}
	@ids = uniq(@ids);
	## @ids
	$ids_for{$cog} = \@ids;
	## %ids_for
	
	return %ids_for;
}

my $end = time; 

print "Total time : ", $end - $begin;
 	

