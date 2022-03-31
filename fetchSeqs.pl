#!/usr/bin/perl
# This script was written by Aymeric Naome and Sinaeda Anderssen

use Modern::Perl '2010';
use strict;
use File::Slurp 'slurp';
use Path::Class 'file';
use Bio::SeqIO;
use Bio::DB::Fasta;
use Smart::Comments;
#use Bio::Tools::Run::Alignment::MAFFT;
#use Bio::PrimarySeqI;

my $fastaFileDir = shift;
my $query_files = shift; #COG_SCO0275.fasta
my $fastaoutdir = shift;
my $fastaoutsuffix = shift;
my $mode = shift;
my @coords=@ARGV;


# Checks if there are coordinates as input, otherwise it takes the CDS for alignment (why?)
# this works
if (@coords) {
	## coords:  @coords;
	if (scalar(@coords) != 2) {
		die "Input TWO bounds, ".scalar(@coords)." given";
	}
	if ( ($coords[0] !~ /[+-]\d+$/) && ($coords[0] !~ /start|end/)) {#&& ($coords[0] !~ m/end/)) {
		die "Upstream bound must be either a +/- integer or \"start\" or \"end\" or a combination of both";
	}
        if ( ($coords[1] !~ /[+-]\d+$/)  &&  ($coords[1] !~ /start|end/)) { #&& ($coords[1] !~ m/end/)) {
                die "Downstream bound must be either a +/- integer or \"start\" or \"end\" or a combination of both";
        }
	#say "$coords[0] $coords[1]";
	my @coords_1 = $coords[0] =~  /([a-z]*)([+-]?\d*)/; # split string and numeric parts 
	my @coords_2 = $coords[1] =~  /([a-z]*)([+-]?\d*)/;
	@coords_1 = map { $_ eq '' ? $_ = "0" : $_ } @coords_1; #replace empty strings with zero 
        @coords_2 = map { $_ eq '' ? $_ = "0" : $_ } @coords_2;
	if ($coords_1[0] ne 0) { if ($coords_1[0] eq "start") { $coords_1[0]=0} else {$coords_1[0]=0.5} } # set the value of "start" position to "0" and "end" position to "0.5"
        if ($coords_2[0] ne 0) { if ($coords_2[0] eq "start") { $coords_2[0]=0} else {$coords_2[0]=0.5} }
	$coords[0] = 0; map {$coords[0]+=$_} @coords_1; # sum the start/end and numerical values
        $coords[1] = 0; map {$coords[1]+=$_} @coords_2;
	#say "$coords[0] $coords[1]";
        ## 1+2: @coords
        say $coords[0]-0.5;
	#say $coords[1]-0.5;
	if ($coords[1] <= $coords[0]) {
		die "Downstream bound must be greater than upstream bound";
	}
	say "Will retrieve sequences from position $coords[0] to position $coords[1] of CDSs for alignment.";
}
else {
	say "No bounds given, will retrieve CDSs for alignment.";
}

# fasta2hash takes the cds.faa files, and creates a hash containing all sequences and their headers, referenced to the seqid number
# Returns a hash
sub fasta2hash {
        my $file = shift;
        my $cut = shift;
        if (!defined($cut)) {$cut = 0;}
        my %hash;
        open(FILE,"<$file") || die("fasta2hash(): Could not open file: $file\n");
        my $last_hash = "";
        while(<FILE>) {
                $_ =~ s/\r/\n/g;
                chomp;
                if ($_ =~ /^>/) {
                        $_ =~ s/^>//;
                        my $header = $_;
                        $_ =~ s/ .*$//s;
                        if ($cut) {
                                $header = $_;
                        }
                        $hash{$_} = [$header,""];
                        $last_hash = $_;
                        next;
                }
                $hash{$last_hash}[1] .= $_."\n";
        }
        close(FILE);
        return %hash;
}
sub modbounds {
        my $start = shift;
        my $end = shift;
        my $strand = shift;
        my $start2;
        my $end2;
if ($strand eq "-") { ($start,$end) = ($end,$start) }
        if (@coords) {
	if ( ($coords[0]-0.5) =~ /^-?\d+$/) {$start2=$end+($strand."1")*($coords[0]-0.5)}
	else {$start2=$start+($strand."1")*$coords[0]}
                if ( ($coords[1]-0.5) =~ /^-?\d+$/) {$end2=$end+($strand."1")*($coords[1]-0.5)}
                else {$end2=$start+($strand."1")*$coords[1]}
        }
        else {$start2=$start; $end2=$end}

        return ($start2,$end2)
}
sub strand_arrow {
	my $strand = shift;
	my $face = shift;
	if (!defined($face)) {$face="-"}
	my $arrow;
	if ($strand eq "-") {$arrow = "<$face$face$face"}
	if ($strand eq "+") {$arrow = "$face$face$face>"}
	if ($strand eq "x") {$arrow = "||"}
	return $arrow
}

# Creating the database of fna files
# if it's the first time it runs on a new set of .fna files, I should de-comment -clean and reindex, and then re-comment them after.
my $db = Bio::DB::Fasta->new($fastaFileDir, -glob => "*.fna", -debug => 1);#, -clean => 1, reindex => 1);

# Creating hash %proteins_for with fasta2hash subroutine
my @protein_files = glob "$fastaFileDir/*.cds.faa";
my %proteins_for;

## The problem was this regex: my ($species) = $protein_file =~ m/(\d+)/;
# it was too "general", only taking a number.
# but for me: $protein_file was "~/Documents/bioinfo/proteinortho180/0.cds.faa"
# this regex thus took the first number it could find, as regexes are usually lazy
# hence: it took 180, and kept that as the single key. Thus, overwriting each value for the corresponding key in the hash, each time.
# So in the end, our hash only contained the key 180, and the values were the hash for the last genome that was processed, which was 99.cds.faa
for my $protein_file (@protein_files) {
    my ($species) = $protein_file =~ m/(\d+).cds.faa/;
	my %hash=fasta2hash($protein_file,0);
	$proteins_for{$species} = \%hash;
    my @keys = keys %proteins_for;
}

open my $justify, '>', "justify.out";
my @queryFiles;
if (! -f $query_files) {
	@queryFiles = glob "$query_files";
	die ("No such file(s), $query_files") if scalar(@queryFiles) == 0
}
else {
	say "Trying to parse file \"$query_files\" for file names";
	@queryFiles = file($query_files)->slurp;
	chomp $queryFiles[0];
	if (! -f $queryFiles[0]) {
		say "File \"$query_files\" does not contain file names, treat as fasta file(s)";
		@queryFiles = glob "$query_files";
	}
}

for my $queryFile (@queryFiles) {
say "Working on file $queryFile ...";
my $basename = (split /\./ ,(split /\//, $queryFile)[-1])[0];
### $basename
mkdir $fastaoutdir;
my $fastaout=$fastaoutdir."/".$basename."_".$fastaoutsuffix.".fasta";
my $out = Bio::SeqIO->new(-file => ">$fastaout" , '-format' => 'FASTA');

open (IN, $queryFile);
LINE:
while (<IN>){
	chomp;
	next if !/^\>/;
	my @fields = split(" ", $_);
	my @location = grep(/location/i, @fields);
	my ($start,$end) = $location[-1] =~ m/[><]?(\d+)[.><]+(\d+)/g;
    my $strand =  $location[0] =~ /complement/ ? "-" : "+";
	my ($seq) = $fields[0] =~ m/^\>([\d_]+)/;
	#my $seqid_aa = $db->get_Seq_by_id($seq);
	#my @ids      = $db->get_all_primary_ids;
	#say @ids;i
	my $species = (split "_", $seq)[0];
	my ($species_code) = $fields[1] =~ /lcl\|(.*)_cds/;
    my $seqid_nt = $db->get_Seq_by_id($species_code);

#	my $seqid_nt = $db->get_Seq_by_id($species);
	my $ntlength =  $seqid_nt->length();
	#say "sequence length " . $ntlength;
	if (!defined( $seqid_nt )) {
                say "Sequence $species not found.";
                next LINE;
        }

    my $seqid_aa = (split " ",$proteins_for{$species}{$seq}[0])[0];
    my @seqid_test = (split " ",$proteins_for{$species}{$seq}[0]);

    if (!defined( $seqid_aa )) {
		say "Sequence $seq not found.";
        next LINE;
    }
	if ($mode eq "prot") {
		$out->write_seq($seqid_aa);
		next LINE;
	}
	my $prevseq = join "_", $species, (split "_", $seqid_aa)[1]-1;
	my $nextseq = join "_", $species, (split "_", $seqid_aa)[1]+1;
	#say $proteins_for{$species}{$prevseq}[0];
	my ($location,$prevend) = $proteins_for{$species}{$prevseq}[0] =~ m/.*location(.*)[=\(><]\d+[.><]+(\d+)/m;
	my $prevstrand;
	if (!defined($prevend)) {
		#say "No upstream CDS";
		$prevend = 0;
		$prevstrand="x";
		
	}
	else {$prevstrand = $location  =~ /complement/ ? "-" : "+";}
	#say $proteins_for{$species}{$nextseq}[0];
	my ($location,$nextstart) = $proteins_for{$species}{$nextseq}[0]  =~ m/.*location(.*)[=\(><](\d+)[.><]+\d+/m;
	my $nextstrand;
        if (!defined($nextstart)) {
                #say "No downstream CDS";
                $nextstart = $ntlength;
		$nextstrand="x";
        }
	else {$nextstrand = $location  =~ /complement/ ? "-" : "+";}
	#say $seq.": strand ".$strand."\n(prev.) ...".strand_arrow($prevstrand).$prevend."___".$start.strand_arrow($strand,"=").$end."___".$nextstart.strand_arrow($nextstrand)."... (next)";
	if ($prevend > $start) { #(($prevend > $start) && ($prevstrand eq $strand)) {
		say "warning: previous CDS overlap";
		#next LINE if $mode eq "CDS";
	}
	if ($nextstart < $end) { # (($nextstart < $end) && ($nextstrand eq $strand )) {
		say "warning: next CDS overlap";
                #next LINE if $mode eq "CDS";
	}
	#say "start: $start, end: $end";
	my @modbounds = modbounds($start,$end,$strand);
	my $start2 = shift @modbounds;
	my $end2 = shift @modbounds;
	#say "$start--$end:$strand";
	#say "$start2--$end2\n";

	if ($mode eq "CDS") {
		if ($strand eq "-") {
        	       	if ( ($start2 > $nextstart) ) {$start2=$nextstart-1}
                	if ($end2 < $prevend) {$end2=$prevend+1}
			say "TEST1" if ($start2 < $end2);
			next LINE if ($start2 < $end2) 
		}
		else {
                	if ($start2 < $prevend) {$start2=$prevend+1}
	                if ($end2 > $nextstart) {$end2=$nextstart-1}
                        say "start2 : $start2 & end2 : $end2";
                        say "TEST2" if ($start2 > $end2);
                        next LINE if ($start2 > $end2)
		}
	}
	# For all the sequences in which it doesn't work, it has a weird start coordinate, and an end coordinate of "1" ! WHY?
	# it also always says strand: + while sometimes this is not true
	
	#say "neighbors: prev. end ".strand_arrow($prevstrand).$prevend.", next start ".$nextstart.strand_arrow($nextstrand);

	#say "from ".$start2." to ".$end2;
	my $subseq;
	if ($strand eq "-") {
		$subseq = $db->seq($seqid_nt, $start2 => $end2);
	}
	else {
		$subseq = $seqid_nt->subseq($start2 => $end2);
	}
	if (length($subseq) < 12) { say "warning: length: ".length($subseq); next LINE};
        say "length: ".length($subseq);
	my $seqout = Bio::Seq->new( -display_id => $seqid_aa."(".$coords[0].",".$coords[1].",".length($subseq).")",
                             -seq => $subseq);
	$out->write_seq($seqout);
	say {$justify} sprintf '%500s', $subseq;
	say $subseq;# "$revsubseq\n";
}
}
