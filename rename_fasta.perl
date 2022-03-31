#!/usr/bin/env perl
use warnings;
use Modern::Perl '2011';
use strict;
use File::Slurp 'slurp';
use Path::Class 'file';
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::Fasta;
#use Smart::Comments;

my $organism = shift;
my $seqinfile = shift;
my $headerinfile = shift;
my $gbkfile = shift;
my $outfile = shift;

#my $fasta = Bio::SeqIO->new( '-format' => 'Fasta' , -file => $seqinfile);
my $fastadb = Bio::DB::Fasta->new( $seqinfile );
### $fastadb
my $header = Bio::SeqIO->new( '-format' => 'Fasta' , -file => $headerinfile);

my $out = Bio::SeqIO->new( -format => 'Fasta', -file => ">$outfile");
$gbkfile = file($gbkfile)->slurp;

say $organism;
my @fasta_ids = $fastadb->get_all_primary_ids;
my $i=0;
while ( my $cds = $header->next_seq() ) {
	if ( $cds->description =~ /protein_id/) {
		my ($WP) = $cds->description =~ /protein_id=([.\w]+)/;
		### $WP
		my ($locus_tag) = $cds->description =~ /locus_tag=([.\w]+)/;
		### $locus_tag
		my ($old_locus_tag) = $gbkfile =~ /locus_tag=\"$locus_tag\" \n \s+ \/old_locus_tag=\"(.\w+)/xm;
		### $old_locus_tag
		# my $old_locus_tag = `grep -A1 $locus_tag $gbkfile | head -n2 | grep -E -o "old_locus_tag=\"\w+"` | cut -d"\"" -f2`;
		# say $locus_tag." # ".$old_locus_tag;
		if (!defined $old_locus_tag) {
			$old_locus_tag = "none";
		}
		#say $locus_tag." # ".$old_locus_tag;
		my $fasta_seq = $fastadb->get_Seq_by_id($WP);
		### $fasta_seq
		my $mod_fasta = Bio::Seq->new( -display_id => $organism."_".$i, -description => $cds->id." ".$cds->description." [old_locus_tag=".$old_locus_tag."]",
                             -seq => $fasta_seq->seq());
		#say $seqstr;
		$out->write_seq($mod_fasta);
		$i++
		
	}
}
		
		
		
