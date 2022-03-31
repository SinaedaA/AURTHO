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


# this is a working example of how to call the fetchSeqsAli script, for debugging purposes if I don't remember how it works
# query_files are the fna files, fetchSeqs makes a db with the fasta files
# so the directory should contain all the fna files

### perl pathTo/Script pathTo/num.fna pathTo/OG.fasta pathTo/outputDir [options]
perl ~/Documents/bioinfo/perl_scripts/fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS noCDS300start50 noCDS -300 start+50
perl ~/Documents/bioinfo/perl_scripts/fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS CDS300start50 CDS -300 start+50
perl ~/Documents/bioinfo/perl_scripts/fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS noCDS300 noCDS -300 -1
perl ~/Documents/bioinfo/perl_scripts/fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS CDS300 CDS -300 -1
perl ~/Documents/bioinfo/perl_scripts/fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS noCDS500 noCDS -500 -1
perl ~/Documents/bioinfo/perl_scripts/fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS CDS500 CDS -500 -1
perl ~/Documents/bioinfo/perl_scripts/fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS noCDS100start50 noCDS -100 start+50
perl ~/Documents/bioinfo/perl_scripts/fetchSeqs.pl $path "$2".fasta $pathToTF/"$2"_UPS CDS100start50 CDS -100 start+50


cd $pathToTF/"$2"_UPS/
echo $PWD

for l in `cat ~/Documents/bioinfo/shell_scripts/samples.ids`; do tpage --define basename=$l --define species=$2 ~/Documents/bioinfo/shell_scripts/template_meme.tt > ./$l.meme.sh; done
#for l in `cat ~/Documents/bioinfo/shell_scripts/sample300.ids`; do tpage --define basename=$l --define species=$2 ~/Documents/bioinfo/shell_scripts/template_meme.tt > ./$l.meme.sh; done

for f in *.meme.sh; do chmod a+x $f; done
for f in *.meme.sh; do ./$f; done
