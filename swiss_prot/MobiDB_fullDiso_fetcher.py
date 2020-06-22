#!/usr/bin/env python

import urllib3
import mmap
from tqdm import tqdm

#uniprot_list_file = "HS_UniprotID_reviewed.list"
uniprot_list_file = "xenopus_UniprotID_oneProtPerGene.list"
output_json_file = "xenopus_mobiDB_full.mjson"


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def main():
	with open(output_json_file,"w+") as output:
		with open(uniprot_list_file,"r+") as SwissProtList:
			for UniprotID in tqdm(SwissProtList, total=get_num_lines(uniprot_list_file)):
				url = "https://mobidb.org/ws/%s/disorder" % (UniprotID.rstrip("\n"))
				http = urllib3.PoolManager()
				request = http.request('GET',url)
				output.write(request.data.decode("utf-8")+"\n")


if __name__ == "__main__":
	main()


