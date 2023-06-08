#!/usr/bin/env python3
import sys
from ete3 import NCBITaxa

ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

InputErrMessage = ', Input Error: input only takes a list of integers.\n\n\tUsage example: ./taxid.py 7070 7460 7165\n'
KeyErrMessage = ', NA, Key Error: taxid not in NCBITaxa'

taxids = []
for taxid in sys.argv[1:]:
    try:
        int(taxid)
        taxids.append(taxid)
    except Exception: print(str(taxid)+InputErrMessage)

try:
    taxnames = ncbi.get_taxid_translator(taxids)
except Exception: print(str(taxid)+InputErrMessage)

for taxid in taxids:
    try: print(taxid+',', taxnames[int(taxid)].replace(' ', '_')) ## Additionally replacing any space in species names with an underscore, as some bioinfo tools prefer it.
    except KeyError: print(str(taxid)+KeyErrMessage)
