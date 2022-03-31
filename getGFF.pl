#!/usr/bin/perl

#use Modern::Perl '2011';
use strict;

#my $infile ="/media/vol2/home/anaome/StreptoStrains/proteins/Results_Feb18/WorkingDirectory/SequenceIDs.txt";
#my $infile="excerpt.dat";

#open my $in, '<', $infile;
#open my $out, '>', $infile . "_out";


#my $locustag = "";
#my $genesynonyms = "";
#my $name = "";
#print "*** " . $infile."\n";


my $species=shift;		# index
my $dir=shift;			# directory name (GCF_(...))
#my $source=`grep -P "^$species:" /media/vol2/home/anaome/StreptoStrains/proteins/Results_Feb18/WorkingDirectory/SpeciesIDs.txt`;
#my $source="$dir/$species.cds.faa";

my $infile = "$dir/$species.cds.faa";						# use previously generated X.cds.faa file
open my $in, '<', $infile;
open my $out, '>', $dir."/".$species.".cds.gff";
my $seqcount=1;
while (my $line = <$in>) {
	#print substr($line,0,1);
	next if substr($line,0,1) ne ">";					# go to next line if first character is not >
	chomp $line;
        my @fields = split(" ", $line);						# split into fields according to blank space
	my ($source) = $fields[1] =~ /lcl\|(.*)_cds/;				# take in field 1 (after index), the group between "lcl|" and "_cds"
	#my $prev_species=$species;
	#$species=(split "_", $fields[0])[0];
	#print $fields[0]."\n";
        #print substr($fields[0],1)."\n";
        pop @fields;								# returns last element of array
	my $strand =  $fields[-1] =~ /complement/ ? "-" : "+";			# if in the next to last field, does it find the pattern "complement" ?=yes so strand = "-"
	#print $strand."\n";
	#if ($species!=$prev_species) {
		#open $out, '>', $species . ".gff";
		#$seqcount=1;
		#$source=`grep -P "^$species:" /media/vol2/home/anaome/StreptoStrains/proteins/Results_Feb18/WorkingDirectory/SpeciesIDs.txt`;
	#}
	
	print {$out} join("\t", $source, "ortho", "CDS", $seqcount, $seqcount, ".", $strand, ".", "ID=".substr($fields[0], 1)).";\n"; # print the table with the features
	$seqcount++;
}
