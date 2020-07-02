#!/usr/bin/env python

import urllib3
import mmap
from tqdm import tqdm

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

uniprot_list_file = "humanCdktargets.txt"
output_folder = "humanCDK1targets_results/"



def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def main():
    with open("%sALL_pdbCoverage.tab" % (output_folder, ), "w+") as output:    
        with open(uniprot_list_file, "r+") as SwissProtList:
            for UniprotID in tqdm(SwissProtList, total=get_num_lines(uniprot_list_file)):
                url = "https://www.ebi.ac.uk/pdbe/widgets/unipdb?tsv=1&uniprot=%s" % (UniprotID.rstrip("\n"))
                http = urllib3.PoolManager()
                request = http.request('GET', url)
                output.write(request.data.decode("utf-8") + "\n")


if __name__ == "__main__":
    main()


