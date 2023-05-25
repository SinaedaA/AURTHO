from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffutils
import pandas as pd
import re
import argparse
import os
import glob
import gzip

def main():
    genome_dir, infofile, logfile = parsing_arguments()
    with open(f"{genome_dir}/{logfile}", 'w') as log:
        log.write(f"""### This log file contains the countsof: 
# CDSs (cds_from_genomic), 
# Proteins (protein.faa), 
# CDSs that are proteins (cds.faa created with this program),
# Paralog pairs of proteins (CDSs that are proteins - Proteins)
GCF\tCDS count\tProteins\tCDS-WP\tParalogPairs
""")
    os.chdir(genome_dir)
    curr_wd = os.getcwd()
    info = pd.read_csv(infofile, sep = "\t", header = None)
    for i in range(len(info)):
        strain, gcf = info.loc[i, 0:1]
        os.chdir(gcf)
        n_cds, n_prot, n_cds_prot = write_cds_faa(strain)
        with open(f"../{logfile}", 'a') as log:
            log.write(f"{gcf}\t{n_cds}\t{n_prot}\t{n_cds_prot}\t{n_cds_prot - n_prot}\n")
        write_cds_gff(strain, gcf)
        os.chdir(curr_wd)

def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='reorder_fasta',
        description='Re-orders the protein.faa record for downloaded genomes, according to the order found in cds_from_genomic.fna')

    parser.add_argument('--genome_dir', default = "./genomes",
                        help="Path to the directory containing assembly directories (GCF_somethingsomething/).")
    parser.add_argument('--infofile', default = "./assembly.info",
                        help="Path to assembly.info file that was created when downloading genomes")
    parser.add_argument('--logfile', default = "./logfile", 
                        help="Log file to show if index.cds.faa has been created")

    args = vars(parser.parse_args())
    genome_dir = args["genome_dir"]
    infofile = args["infofile"]
    logfile = args["logfile"]

    return genome_dir, infofile, logfile

def fasta_gz_to_dict(pattern):
    """
    This function takes a filehandle of a fasta file, then reads the content and stores it in a dictionary, which it returns.
    It can handle gz zipped files, as well as unzipped fasta files, and doesn't care about the used alphabet.
    """
    fh = glob.glob(pattern)[0]
    if fh.endswith(".gz"):
        with gzip.open(fh, mode='rt') as zipf:
            content = SeqIO.to_dict(SeqIO.parse(zipf, "fasta"))
    else:
        with open(fh, mode = 'rt') as f:
            content = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    return(content)

def write_cds_faa(strain):
    """
    Takes the assembly.info file as a data.frame, and goes into each assembly directory (GCF something-something) to create a cds.faa file.
    The strain.cds.faa file is made by comparing the cds_from_genomic.fna and the protein.faa file.
    cds_from_genomic.fna = genes that form coding sequences, ATGC alphabet (fna)
    protein.faa = accessioned protein products annotated on the genome assembly (amino acid alphabet)
    strain.cds.faa = headers from cds_from_genomic, with protein sequences from protein.faa
    This file is not the same as translated.faa, as this file is merely a naive translation of cds_from_genomic.
    This function, from the genomes/ directory, goes into each assembly dir and reads the contents of the .gz files of interest, does what it needs, then goes back in /genomes.

    :param info_df: pandas dataframe from assembly.info file
    :param logfile: filehandle for the logfile, which will contain the numbers of CDSs, proteins, and CDSs that are proteins.
    """
    proteins = fasta_gz_to_dict(pattern = f"*protein.faa*")
    cds = fasta_gz_to_dict(pattern = f"*cds_from_genomic.fna*")
    cds_faa = {}
    for key, value in cds.items():
        ## this works for fungi at least, should check if it still works for bacteria from RefSeq with the WP_ notation
        match = re.search(
            "\[protein_id=([A-Z]*_?[0-9]*\.[0-9]{1})\]", value.description
        )
        # match = re.search(
        #     "\[protein_id=(WP_[0-9]*\.[0-9]{1})\]", value.description)
        if match is not None:
            pi = match.group(1)
            cds_faa[key] = SeqRecord(
                seq=proteins[pi].seq, id=f"{pi} {value.id}", name=value.name, description=value.description)
        else:
            continue
    pi_from_cds_faa = re.findall(
        "cds_([A-Z]*_?[0-9]*\.[0-9]{1})", str(cds_faa.keys())
    )
    # wp_from_cds_faa = re.findall(
    #     "(WP_[0-9]*\.[0-9]{1})", str(cds_faa.keys()))
    n_cds = len(cds.keys())
    n_prot = len(proteins.keys()) # same as unique items in wp_from_cds_faa, meaning the additional sequences in latter are paralog proteins
    n_cds_prot = len(pi_from_cds_faa)
    with open(f"{strain}.cds.faa", 'w') as fasta:
        SeqIO.write(cds_faa.values(), fasta, "fasta")
    return n_cds, n_prot, n_cds_prot

def write_cds_gff(strain, gcf):
    gff = glob.glob(f"*genomic.gff*")[0]
    db = gffutils.create_db(data = gff,
                            dbfn=":memory:",
                            merge_strategy='create_unique',
                            disable_infer_transcripts=True,
                            disable_infer_genes=True)
    ## write the *.cds.gff file
    with open(f"{strain}.cds.gff", 'w') as gff_fh:
        for feat in db.all_features(featuretype='CDS'):
            if not re.search("pseudo=true", str(feat)):
                gff_fh.write(f"{str(feat)}\n") 

if __name__ == '__main__':
    main()