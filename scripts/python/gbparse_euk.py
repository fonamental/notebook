#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Blast import NCBIXML
#Usage: $python gbparse_euk.py outfile.txt infile.gb
OUT = open(sys.argv[1], 'w')
OUT.write("Accno\tLength\tOrganism\tTaxonomy\tStrain\tIsolate\tCulture_collection\tVoucher\tIsolation_source\tClone\tOrganelle\tGene\tProduct\tMolecule_type\tNote\tDB_xref\tPrimers\tCountry\tYear_submitted\tAuthor\tKeyword\tSequence\n")
result_handle = open(sys.argv[2])
gbfiles = SeqIO.parse(result_handle, 'gb')
for rec in gbfiles:
    acc = rec.id
    sequence = str(rec.seq)
    length = str(len(rec.seq))
    recfeat1 = rec.features[1]
    source = rec.features[0]
    if 'organism' in rec.annotations:
        organism = rec.annotations['organism']
    if 'taxonomy' in rec.annotations:	
        taxonomy = "_".join(rec.annotations['taxonomy'])
    if 'clone' in source.qualifiers:
        clone = source.qualifiers['clone'][0]
    else:
        clone = ""
    if 'source' in rec.annotations:
        gensource = rec.annotations['source']
    else:
        gensource = ""
    if 'isolation_source' in source.qualifiers:
        isolation_source = source.qualifiers['isolation_source'][0]
    else:
        isolation_source = ""
    if 'country' in source.qualifiers:
        country = source.qualifiers['country'][0]
    else:
        country = ""
    if 'strain' in source.qualifiers:
        strain = source.qualifiers['strain'][0]
    else:
        strain = ""
    if 'isolate' in source.qualifiers:
        isolate = source.qualifiers['isolate'][0]
    else:
        isolate = ""
    if 'culture_collection' in source.qualifiers:
        cultcol = source.qualifiers['culture_collection'][0]
    else:
        cultcol = ""
    if 'specimen_voucher' in source.qualifiers:
        voucher = source.qualifiers['specimen_voucher'][0]
    else:
        voucher = ""
    if 'PCR_primers' in source.qualifiers:
        primers = ".".join(source.qualifiers['PCR_primers'])
    else:
        primers = ""
    if 'gene' in recfeat1.qualifiers:
        whichgene = recfeat1.qualifiers['gene'][0]
    else:
        whichgene = ""
    if 'product' in recfeat1.qualifiers:
        product = recfeat1.qualifiers['product'][0]
    else:
        product=""
    if 'mol_type' in source.qualifiers:
        molecule_type = source.qualifiers['mol_type'][0]
    else:
        molecule_type = ""
    if 'organelle' in source.qualifiers:
        organelle = source.qualifiers['organelle'][0]
    else:
        organelle = ""
    if 'note' in source.qualifiers:
        note = source.qualifiers['note'][0]
    else:
        note = ""
    if 'db_xref' in source.qualifiers:
        db_xref=source.qualifiers['db_xref'][0]
    else:
        db_xref=""
    if 'references' in rec.annotations:
        pubref = rec.annotations['references'][0]
        authors = pubref.authors
        firstaut = authors.split(".,")[0]
    else:
        firstaut = ""
    if 'date' in rec.annotations:
        date = rec.annotations['date']
        submyear = rec.annotations['date'][7:11]
    if 'keywords' in rec.annotations:
        keyword = rec.annotations['keywords'][0]
    fields = [acc, length, organism, taxonomy, strain, isolate, cultcol, voucher, isolation_source, clone, organelle, whichgene, product, molecule_type, note, db_xref, primers, country, submyear, firstaut, keyword, sequence]
    OUT.write("\t".join(fields)+ "\n")
OUT.close()			
