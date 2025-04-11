import os
import argparse
from pathlib import Path
import gzip
import glob
from datetime import date
import pandas as pd
import re

def main():
    outdir, genus, level, overwrite, db, repo = parsing_arguments()
    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)
    info = Path("./assembly.info")
    info.unlink(missing_ok=True)
    month = date.today().strftime("%b%Y")
    if not os.path.exists(f"assembly_summary_{month}.txt"):
        print(f"###### Downloading latest assembly_summary.txt (bacterial) from NCBI FTP server as assembly_summary_{month}.txt")
        os.system(f"rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/{db}/{repo}/assembly_summary.txt ./assembly_summary_{db}_{month}.txt")
    else:
        print(f"##### assembly_summary (bacterial) for {month} already exists.")
    
    print(f"Reading assembly_summary, and writing ftp paths for selected {genus} and {level}")
    ftp_list = filter_assembly(file = f"assembly_summary_{db}_{month}.txt", a_level = level, genus = genus, date = month)
    download_gcf(ftp_list, overwrite)
    make_info_file(f"assembly_summary_{db}_{month}.txt", ftp_list)

def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='download_genomes.py',
        description='Downloads GCF records from the ncbi ftp server for bacteria')

    parser.add_argument('--outdir', default = "./genomes", help = "Path to output directory. Defaults to current_dir/genomes.")
    parser.add_argument('--genus', help = "Specify bacterial genus of interest", required=True)
    parser.add_argument('--level-of-assembly', dest='level', help = "Level of assembly for download", choices = ['Complete Genome', 'Chromosome', 'Contig', 'Scaffold', 'all'], default = "Complete Genome")
    parser.add_argument('--overwrite', help = "Overwrite assembly files if they already exist", default = "True", choices = ["True", "False"])
    parser.add_argument('--database', help = "Download either from 'refseq' or 'genbank'", default = 'refseq', choices = ['refseq', 'genbank'])
    parser.add_argument('--repository', help = 'Which organisms to download', choices = ['bacteria', 'fungi', 'invertebrate', 'metagenomes', 'other', 'plant', 
                                                                                      'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral'])
    args = vars(parser.parse_args())
    outdir = args["outdir"]
    genus = args["genus"]
    level = args["level"]
    db = args["database"]
    repo = args["repository"]
    overwrite = args["overwrite"]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f"{outdir} did not exist, creating it.")
    return outdir, genus, level, overwrite, db, repo

def filter_assembly(file, a_level, genus, date):
    assembly_file = pd.read_csv(file, sep="\t", header=1, low_memory=False)
    pd.options.display.max_colwidth = 500
    if a_level != "all":
        ftp_path = assembly_file[(assembly_file["assembly_level"] == a_level) & (
            assembly_file["organism_name"].str.contains(genus))]["ftp_path"].to_string(index=False)
    elif a_level == "all":
        ftp_path = assembly_file[(assembly_file["organism_name"].str.contains(genus))]["ftp_path"].to_string(index = False)
    ftp_path = re.sub("https://", "", ftp_path)
    ftp_list = list(map(str.strip, ftp_path.split("\n")))
    with open(f"ftpdirpaths_{date}.txt", mode ='wt', encoding = "utf-8") as myftp:
        myftp.write('\n'.join(ftp_list))
    return ftp_list

def download_gcf(ftp_list, overwrite):
    #for i in range(10):
    for i in range(len(ftp_list)):
        gcf = ftp_list[i].split("/")[-1]
        if overwrite == "True":
            print(f"Downloading assembly {gcf}")
            os.system(
                f"rsync --copy-links --times --verbose -r --keep-dirlinks rsync://{ftp_list[i]} ./")
        else:
            if os.path.exists(gcf):
                print(
                    f"Directory for assembly {gcf} already exists, skipping...")
                continue
            else:
                print(f"Downloading assembly {gcf}")
                os.system(
                    f"rsync --copy-links --times --verbose -r --keep-dirlinks rsync://{ftp_list[i]} ./")

def make_info_file(file, ftp_list):
    #for i in range(10):
    for i in range(len(ftp_list)):
        assembly_file = pd.read_csv(file, sep="\t", header=1, low_memory=False)
        gcf = ftp_list[i].split("/")[-1]
        print(gcf)
        gbff = glob.glob(f"./{gcf}/*gbff*")[0]
        organism = assembly_file[assembly_file["ftp_path"].str.contains(gcf)]["organism_name"].to_string(index = False)
        if gbff.endswith(".gz"):
            with gzip.open(gbff, mode = 'r') as zipf:
                lt, olt = get_lt_olt(zipf, zip = True)
        else:
            with open(gbff, mode = 'r') as f:
                lt, olt = get_lt_olt(f, zip = False)
        with open("assembly.info", 'a') as info:
            info.write(f"{i}\t{gcf}\t{organism}\t{lt}\t{olt}\n")

def get_lt_olt(filehandle, zip = True):
    if zip == True:
        content = filehandle.read().decode("utf-8")
    else:
        content = filehandle.read()
    locustag = re.search(
        "(?i)/locus_tag=\"([A-Z0-9]*)_*", content)
    if locustag is not None:
        lt = locustag.group(1)
    else:
        lt = ""
    old_locustag = re.search(
        "(?i)/old_locus_tag=\"([A-Z0-9]*)_*", content)
    if old_locustag is not None:
        olt = old_locustag.group(1)
    else:
        olt = ""
    return(lt, olt)

if __name__ == '__main__':
    main()
