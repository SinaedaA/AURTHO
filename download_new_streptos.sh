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

for dir in `cat ftpdirpaths_$date`; do outdir=`echo $dir | rev | cut -d"/" -f1 | rev`; ftp=`echo $dir | cut -d"/" -f3-100`; assembly=`echo $outdir | cut -d"_" -f1-2`; ; done

ftp=`echo $dir | cut -d"/" -f3-100`; assembly=`echo $outdir | cut -d"_" -f1-2`; 
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