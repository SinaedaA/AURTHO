# AURTHO methodology
----------------------------------------------

## 1. Download of latest _Streptomyces_ genomes

```
cd ~/Documents/bioinfo/strepto-assembly

../shell_scripts/download_new_streptos.sh
```
This script will look inside the strepto-assembly directory, and only download genomes with GCF ids that are not yet present inside this directory. 
It will create a new directory named `strepto_genomes_$date`, as well as 2 files: `assembly_summ_$date` and `ftpdirpaths_$date` (containing the links used to download the genomes).  
Here is the [link](#download_new_streptos) to the script that is used for the download. 
  
### Unzip the files of interest
```
cd strepto_genomes_$date
gunzip */*protein.faa.gz
gunzip */*cds_from_genomic.fna.gz
gunzip */*_genomic.gbff.gz
```
The files that are most interesting for the orthogroup creation are the annotated protein products (.faa), the CDSs (.cds\_from\_genomic.fna), and the GenBank flat file annotated genome sequence (.gbff, what we see on NCBI).  

### Creating the assembly.info file
This file contains all the information on each strain, including the locus\_tag, the old\_locus\_tag, the index that we gave the strain (which is given according to the order in which it appears inside the directory).  
The **i** should be set to 0 (if it is the first time we download genomes), or to _the number of strains last downloaded_.  

I recently have had to change the printf statement to an echo, in order to get each line on a newline. The "\t" is replaced by "<ctrl v+i>" in the command line.

```
rm −f assembly.infoi=0

for dir in `ls -d GCF_* | cut -d"/" -f1`;	do organism=`awk -F "\t" '$1=="'$(echo $dir | cut -d"_" -f1-2)'"{print $8,$9,$10}' ../assembly_summ_030718.txt`	locus_tag=`egrep -o -m1 "locus_tag=\"\w+" $dir/*_genomic.gbff | cut -d"\"" -f2 | cut -d"_" -f1`	old_locus_tag=`egrep -o -m1 "old_locus_tag=\"\w+" $dir/*_genomic.gbff |cut -d"\"" -f2 | cut -d"_" -f1`	printf $i"\t"$dir"\t"$organism"\t"$locus_tag"\t"$old_locus_tag >> assembly.info	i=$(($i+1)) 
done
```  
### Reordering the sequences in \*protein.faa according to the order of the sequences in \*cds\_from\_genomic.fna
```
IFS=$'\n'
for strain in `cat assembly.info`
 do i=`echo $strain | awk '{print $1}'`
 dir=`echo $strain | awk '{print $2}'`
 perl rename_fasta.perl $i $dir/"$dir"_protein.faa $dir/"$dir"_cds_from_genomic.fna $dir/"$dir"_genomic.gbff $dir/$i.cds.faa
 perl getGFF.pl $i $dir
 ln -s $PWD/$dir/$i.cds.faa ~/Documents/bioinfo/proteinortho/$i.cds.faa
 ln -s $PWD/$dir/$i.cds.gff ~/Documents/bioinfo/proteinortho/$i.cds.gff
done

unset IFS
```
### ProteinOrtho

```
proteinortho6.pl -project=streptomyces -singles -cpus=1 -selfblast -verbose ~/path/to/*cds.faa
```
Once this step is done, we copy the *faa files inside the proteinortho/ directory, which will be useful for the grab_protein script.

I modified the grab_protein script, in order to rename all the orthogroups according to the best studied strains. For example, if an orthogroup contains a protein from *Streptomyces coelicolor*, then this orthogroup will be named after the first SCO protein in the group.   

```
mkdir COG

perl grab_protein_sina.pl streptomyces.proteinortho ; find . -name "*fasta" -print0 | xargs -0 -J % mv % COG/
```   
	    
## 2. Identifying LacI transcription factor COGs
The next step is to select only the orthogroups that contain LacI-family transcription factors.

These code chunk will perform an hmmscan on all the *.cds.faa files. Hence, it can be performed at the same time as the clustering step of ProteinOrtho. 
`TFscript.sh` is the script incorporating all the information from the P2TF website, which defines specific domain combinations for each TF family. 
```
for strain in 7 3 2 0 1 13 $(seq 4 6)  $(seq 8 12) $(seq 14 89)
 do
 ~/hmmer3/binaries/hmmscan --cpu 2  --tblout TF_scan_"$strain".out  --domtblout TF_scan_dom_"$strain".out --pfamtblout TF_scan_pfam_"$strain".out ~/Documents/bioinfo/hmm/TF.hmm ../"$strain".cds.faa
 bash ~/Documents/bioinfo/shell_scripts/TFscript.sh TF_scan_"$strain".out > TF_filtered_"$strain".dat
done
```   
#####Extracting sequence IDs for a certain family of TFs (here LacI):

```
for strain in 7 3 2 0 1 13 $(seq 4 6)  $(seq 8 12) $(seq 14 181)
 do egrep -e "^\d+_\d+ LacI"  TF_filtered_"$strain".dat | awk '{print ">"$1" lcl"}' >> LacI_seqids.dat
done
```   
#####Grep the sequence IDs and the COG fasta file names and store them in LacI_COG.dat   
The [find_OGs\_TF.pl](#find_OGs_TF.pl) script outputs which sequences in which COGs contain at least one TF determined to be a LacI TF based on the hmmscan and TFscript.  

`Usage: perl find_OGs_TF.pl LacI_seqids.dat LacI_COG.dat`      


## 3. Functional coherence check for COGs containing LacI-TF
```
IFS=$'\n'
for cog in `cat LacI_COG.dat | cut -d":" -f1 | sort | uniq -c | sort -g -r -k1`
	do a=`echo $cog | awk '{print $2}'`
	count=`grep -c ">" $a`
	descr=`grep ">" $a | perl -ne 'print "$1\n" if /protein=([\w\s-]+)/' | sort | uniq -c | sort -g -r -k1`
	echo $cog $count $descr
done > LacI_temp
``` 

The LacI_temp file looks like this:   

```
  93 ./OG-SCO0456.fasta 95   93 LacI family transcriptional regulator    2 transcriptional regulator
  91 ./OG-SCO2232.fasta 91   85 LacI family transcriptional regulator    3 LacI family DNA-binding transcriptional regulator    1 transcriptional regulator    1 maltose operon transcriptional repressor    1 MalR repressor protein
  90 ./OG-SCO4158.fasta 90   89 LacI family transcriptional regulator    1 LacI-family regulatory protein
  90 ./OG-SCO3943.fasta 90   50 transcriptional regulator   40 LacI family transcriptional regulator
  88 ./OG-SCO1078.fasta 88   88 LacI family transcriptional regulator
```
The first column indicates the number of proteins that were detected by `hmmscan` to have a typical LacI domain.
The thirst column is the total number of sequences in the COG.
The following columns counts the number of times the adjacent annotation was found in the proteins of that COG.

Below are some examples of COGs for which the annotation did not pass the manual check of functional coherence with a regulatory role. 

```
  11 ./OG-SCAB_RS07735.fasta 12    8 N-6 DNA methylase    4 SAM-dependent methyltransferase
   8 ./OG-SGR_RS04355.fasta 50   23 hypothetical protein   11 transcriptional regulator    4 streptomycin biosynthesis regulator    4 streptomycin biosynthesis operon regulator    3 streptomycin biosynthesis protein    3 nuclease    2 Streptomycin biosynthesis operon regulatory protein
   5 ./OG-SACTE_RS28705.fasta 12    7 N-6 DNA methylase    3 hypothetical protein    1 type II restriction endonuclease subunit M    1 restriction endonuclease subunit M
   4 ./OG-SCO2312.fasta 90   67 hypothetical protein   20 class I SAM-dependent methyltransferase    3 secreted 
```   

## 4-5. Extracting upstream regions of LacI-COGs and MEME

Script that extracts upstream regions and performs MEME.
It needs a `samples.ids` file, which contains the parameters for the script that fetches the upstream sequences:

```
CDS100start50
CDS300
CDS300end
CDS300start100
CDS500
```

And a template_meme.tt file, containing the parameters to perform MEME.

```
for file in [% species %]_[% basename %].fasta
	do meme $file -oc [% species %]_meme/[% basename %]_anr_long -dna -mod anr -maxw 30 -nmotifs 4 -revcomp
	meme $file -oc [% species %]_meme/[% basename %]_zoops_long -dna -mod zoops -maxw 30 -nmotifs 4 -revcomp
	meme $file -oc [% species %]_meme/[% basename %]_anr_medium -dna -mod anr -maxw 20 -nmotifs 4 -revcomp
	meme $file -oc [% species %]_meme/[% basename %]_zoops_medium -dna -mod zoops -maxw 20 -nmotifs 4 -revcomp
	meme $file -oc [% species %]_meme/[% basename %]_anr_short -dna -mod anr -maxw 10 -nmotifs 4 -revcomp
	meme $file -oc [% species %]_meme/[% basename %]_zoops_short -dna -mod zoops -maxw 10 -nmotifs 4 -revcomp
done
```

Usage of [UPS\_meme\_auto.sh](#UPS_meme_auto.sh)


`UPS_meme_auto.sh ~/path/to/genomes.fna name_of_COG ~/path/to/output/folder/`

Use it in a `for` loop to run it on all COGs that you want to analyse.

# SCRIPTS
-----------

##<a name="download_new_streptos"></a> Download new _Streptomyces_ genomes script

```
date=`date +"%m%d%y"`

echo "	##### Downloading latest assembly_summary.txt from NCBI FTP server as assembly_summ_"$date".txt"

rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt ./assembly_summ_"$date".txt

awk -F "\t" '$8 ~ /Streptomyces/ && $12=="Complete Genome" && $11=="latest"{print $20}' assembly_summ_"$date".txt > ftpdirpaths_$date

for dir in `cat ftpdirpaths_$date`
	do outdir=`echo $dir | rev | cut -d"/" -f1 | rev`
	if [ ! -d ./"$outdir" ] 
	then
		ftp=`echo $dir | cut -d"/" -f3-100`
		assembly=`echo $outdir | cut -d"_" -f1-2`
		organism=`awk -F "\t" '$1=="'$assembly'"{print $8,$9,$10}' assembly_summ_"$date".txt`
		echo "#### Downloading assembly $assembly ($organism) to ./outdir"
		rsync --copy-links --times --verbose -r --keep-dirlinks rsync://$ftp ./strepto_genomes_"$date"/
	fi
done
```

##<a name="find_OGs_TF.pl"></a> find\_OGs\_TF.pl

```
#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;
use Smart::Comments;
use Tie::IxHash;
use List::AllUtils qw(uniq);
use Storable;

my $begin = time; 

# Open directory and get all the files in it (perl has no 'argument list too long' error)
opendir(DIR, "/Users/sinaeda/Documents/bioinfo/ProteinOrtho/COG90/");
my @infiles = grep(/OG-.*\.fasta/, readdir(DIR));
closedir(DIR);


# Create tf_ids, which opens and then closes the dat file for the TF
# Only once
my $dat = shift;
my $outfile = shift;
open my $dh, '<', $dat;
my @tf_ids;
while (my $l = <$dh>) {
	chomp $l;
	push @tf_ids, $l;
}
close $dh;


# Create hash for COGs and their corresponding gene IDs
tie my %ids_for, 'Tie::IxHash';

if (! -e '/Users/sinaeda/Documents/bioinfo/ProteinOrtho/COG90/COG_hash.ref') {
	for my $infile (@infiles) {
		## $infile
		%ids_for = (%ids_for, read_COG_fasta($infile));
	}
	## %ids_for
	store \%ids_for, '/Users/sinaeda/Documents/bioinfo/ProteinOrtho/COG90/COG_hash.ref';
}

my $id_ref = retrieve('/Users/sinaeda/Documents/bioinfo/ProteinOrtho/COG90/COG_hash.ref');
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
```



##<a name="fetchSeqs.pl"></a> fetchSeqs.pl   
```
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

```


##<a name="UPS_meme_auto"></a> UPS\_meme\_auto.sh

```
### Make script for extracting UPS sequences from different COGs again and again to ease the process

## Get upstream regions
## usage
## for f in `cat group.txt`; do ./test_ups_meme.sh $f group_name; done
## in group.txt should contain name of OGs (e.g. OG-SBI_RS36640) separated by newline

# for f in `cat group.txt`; do group=`echo $f | cut -d'-' -f2`; ./UPS_auto.sh ~/Documents/bioinfo/surfactin/COG/$f $group $f; done

# path = ~/Documents/bioinfo/proteinortho180/
path=$1
OG=$2
pathToTF=$3
#SCO=$3

# Works
mkdir $pathToTF/"$2"_UPS/
mkdir $pathToTF/"$2"_UPS/"$2"_meme

# Go to directory containing ... files
cd $path/COG180

### perl pathTo/Script pathTo/num.fna pathTo/OG.fasta pathTo/outputDir [options]
perl fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS noCDS300start50 noCDS -300 start+50
perl fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS CDS300start50 CDS -300 start+50
perl fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS noCDS300 noCDS -300 -1
perl fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS CDS300 CDS -300 -1
perl fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS noCDS500 noCDS -500 -1
perl fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS CDS500 CDS -500 -1
perl fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS noCDS100start50 noCDS -100 start+50
perl fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS CDS100start50 CDS -100 start+50


cd $pathToTF/"$2"_UPS/
echo $PWD

for l in `cat ~/Documents/bioinfo/shell_scripts/samples.ids`; do tpage --define basename=$l --define species=$2 ~/Documents/bioinfo/shell_scripts/template_meme.tt > ./$l.meme.sh; done
#for l in `cat ~/Documents/bioinfo/shell_scripts/sample300.ids`; do tpage --define basename=$l --define species=$2 ~/Documents/bioinfo/shell_scripts/template_meme.tt > ./$l.meme.sh; done

for f in *.meme.sh; do chmod a+x $f; done
for f in *.meme.sh; do ./$f; done

```