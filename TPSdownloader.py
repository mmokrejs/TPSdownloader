#! /usr/bin/python3

import sys
import os
import datetime
import numpy as np
import pandas as pd
from shlex import split as shlex_split
from subprocess import Popen, PIPE, call #, STDOUT
from optparse import OptionParser
import xml.etree.ElementTree as ET
import time
import re
# from re import compile, findall, finditer, sub, IGNORECASE
import copy
from itertools import chain
import gzip

version = "20210511"
myparser = OptionParser(version="%s version %s" % ('%prog', version))
myparser.add_option("--uniprot-ids-from-file", action="store", type="string", dest="uniprot_ids_from_file", default='',
    help="Obtain a list of Uniprot IDs from column 'Uniprot ID' in 'Sheet1'.")
myparser.add_option("--uniprot-id-file", action="store", type="string", dest="uniprot_idfile", default='-',
    help="Obtain list of Uniprot IDs from a single-column text file or STDIN, is mutually exclusive with --chebi-id-file.")
myparser.add_option("--chebi-id-file", action="store", type="string", dest="chebi_idfile", default='-',
    help="Obtain a list of ChEBI IDs from a single-column text file or STDIN, is mutually exclusive with --uniprot-id-file.")
myparser.add_option("--uniprot-id", action="store", type="string", dest="uniprot_id", default=None,
    help="List of Uniprot IDs, is mutually exclusive with --chebi-id.")
myparser.add_option("--chebi-id", action="store", type="string", dest="chebi_id", default=None,
    help="List of ChEBI IDs, is mutually exclusive with --uniprot-id.")
myparser.add_option("--xls-storage", action="store", type="string", dest="xls_storage", default="TPSdownloader.xls",
    help="Use this file to parse input data and also as an output file with new results appended to it. Use None to disable the default. Default is TPSdownloader.xls.")
myparser.add_option("--outfmt", action="store", type="string", dest="outfmt", default="xls",
    help="Format of output file. It is used to preserve data between restarts too. CSV or XLSX (default)")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug to some value")
(myoptions, myargs) = myparser.parse_args()


extra_product_colnames = ['Name of product', 'Product compound description', 'Chemical formula of product', 'SMILES of product (including stereochemistry)']
extra_substrate_colnames = ['Substrate (including stereochemistry)', 'Substrate compound description', 'Chemical formula of substrate', 'SMILES of substrate (including stereochemistry)']
extra_cofactor_colnames = ['Cofactors', 'Cofactors compound description', 'Chemical formula of cofactor', 'SMILES of cofactor (including stereochemistry)']
extra_intermediate_colnames = ['Name of intermediate', 'Intermediate compound description', 'Chemical formula of intermediate', 'SMILES of intermediate (including stereochemistry)']


def parse_list_or_set_line(line):
    _strlist = line.split(':')[1][1:-1]
    return filter(lambda x: x != 'None' and x is not None, parse_list_or_set(_strlist))


def parse_list_or_set(_strlist):
    if _strlist == 'None':
         _out = []
    elif _strlist == 'set()':
        _out = []
    elif _strlist == '[]':
        _out = []
    elif _strlist and _strlist[0] == '[':
        if "," in _strlist:
            _out = _strlist[1:-1].split(', ')
            if _out[0][0] == "'":
                _out = map(lambda y: y[1:-1], _out)
        elif  _strlist[:2] == "['" or _strlist[:2] == '["':
            _out = [_strlist[2:-2]]
        elif  _strlist == '[None]':
            _out = [None]
        else:
            _out = [_strlist[2:-2]]

    elif _strlist and _strlist[0] == 'set(': # set(['aa', 'a', 'c', 'd', 'bb'])
        _out = map(lambda x: x[1:-1], _strlist[5:-2].split(', '))
    else:
        _out = _strlist

    # convert to integers if possible
    if isinstance(_out, list):
        _new = []
        for _item in _out:
            try:
                _new.append(int(_item))
            except:
                _new.append(_item)

        return _new
    else:
        return _out


def parse_known_terpenes(filename="terpene_names.uniq.txt"):
    """Multi-line strings are wrapped by double-quotes. Such entries appear in the XLSX files
    mostly by mistake but line-based approach yields entries starting with double-quote sign.
    """

    _terpenes = []
    with open("terpene_names.uniq.txt") as _file:
        for _line in _file:
            if _line[0] != '#':
                if _line[0] != '"':
                    _terpenes.append(_line[:-1])
                else:
                    _terpenes.append(_line[1:-1])
    if myoptions.debug: print("Debug: terpenes=%s" % str(_terpenes))
    return _terpenes


# blacklist of IDs of reaction substrates or non-cyclic terpenes, etc., but *also* including β-farnesene CHEBI:10418
# these IDs are not downloaded into the cache, unless they already were downloaded before addition to this list
# IDs appearing in this list also do not get output into the output list of terpenes
non_terpene_chebi_ids = ['CHEBI:35194', 'CHEBI:15377', 'CHEBI:33019', 'CHEBI:58057', 'CHEBI:128769', 'CHEBI:175763', 'CHEBI:57533', 'CHEBI:10418', 'CHEBI:10280', 'CHEBI:175763', 'CHEBI:64283', 'CHEBI:58756', 'CHEBI:15441', 'CHEBI:58622', 'CHEBI:58553', 'CHEBI:57665', 'CHEBI:58635', 'CHEBI:15440', 'CHEBI:138223', 'CHEBI:128769', 'CHEBI:57907', 'CHEBI:64801', 'CHEBI:61984', 'CHEBI:58206', 'CHEBI:138167', 'CHEBI:15347', 'CHEBI:162247', 'CHEBI:17221', 'CHEBI:60374', 'CHEBI:61746', 'CHEBI:98', 'CHEBI:46702', 'CHEBI:61987', 'CHEBI:16240', 'CHEBI:35757', 'CHEBI:3407', 'CHEBI:13657', 'CHEBI:25382', 'CHEBI:43474', 'CHEBI:43470', 'CHEBI:29139', 'CHEBI:61987', 'CHEBI:28938', 'CHEBI:83628', 'CHEBI:24646', 'CHEBI:134188', 'CHEBI:28938', 'CHEBI:22534', 'CHEBI:49783', 'CHEBI:7435', 'CHEBI:139521', 'CHEBI:15379', 'CHEBI:44742', 'CHEBI:7860', 'CHEBI:10745', 'CHEBI:13416', 'CHEBI:23833', 'CHEBI:25366', 'CHEBI:29097', 'CHEBI:30491', 'CHEBI:139520', 'CHEBI:132124', 'CHEBI:57540', 'CHEBI:58340', 'CHEBI:128753', 'CHEBI:33384', 'CHEBI:17268', 'CHEBI:57288', 'CHEBI:33738', 'CHEBI:33737', 'CHEBI:58720', 'CHEBI:57783', 'CHEBI:57287', 'CHEBI:15378', 'CHEBI:57623']
# CHEBI:35194 - isoprene
# CHEBI:33019 - diphosphate(3−)
# CHEBI:58057 - geranyl diphosphate(3−)
# CHEBI:57623 - prenyl diphosphate(3-)
# CHEBI:128769 - isopentenyl diphosphate(3−)
# CHEBI:175763 - 2-trans,6-trans-farnesyl diphosphate(3−)
# CHEBI:57533 - geranylgeranyl diphosphate(3−)
# CHEBI:15377 - water
# CHEBI:16240 - hydrogen peroxide
# CHEBI:35757 - monocarboxylic acid anion
# CHEBI:3407 - monocarboxylic acid anion
# CHEBI:13657 - monocarboxylic acid anion
# CHEBI:25382 - monocarboxylic acid anion
# CHEBI:43474 - hydrogenphosphate
# CHEBI:43470 - hydrogenphosphate
# CHEBI:29139 - hydrogenphosphate
# CHEBI:61987 - 2-methylisoborneol
# CHEBI:28938 - ammonium
# CHEBI:83628 - N-acylammonia
# CHEBI:24646 - hydroquinones
# CHEBI:134188 - hydroquinones
# CHEBI:28938 - ammonium
# CHEBI:22534 - ammonium
# CHEBI:49783 - ammonium
# CHEBI:7435 - ammonium
# CHEBI:139521 - phenolic radical donor
# CHEBI:15379 - dioxygen
# CHEBI:44742 - dioxygen
# CHEBI:7860 - dioxygen
# CHEBI:10745 - dioxygen
# CHEBI:13416 - dioxygen
# CHEBI:23833 - dioxygen
# CHEBI:25366 - dioxygen
# CHEBI:29097 - dioxygen
# CHEBI:30491 - dioxygen
# CHEBI:139520 - phenolic donor
# CHEBI:132124 - 1,4-benzoquinones
# CHEBI:57540 - NAD(1-)
# CHEBI:58340 - O-acetyl-L-serine zwitterion
# CHEBI:128753 - (2E)-4-hydroxy-3-methylbut-2-enyl diphosphate(3-)
# CHEBI:33384 - L-serine zwitterion
# CHEBI:17268 - myo-inositol
# CHEBI:57288 - acetyl-CoA(4-)
# CHEBI:57287 - coenzyme A(4-)
# CHEBI:33738 - di-mu-sulfido-diiron(1+)
# CHEBI:33737 - di-mu-sulfido-diiron(2+)
# CHEBI:58720 - D-glucopyranuronate
# CHEBI:57783 - NADPH(4-)
# CHEBI:15378 - hydron





def parse_chebi_xml(filename):
    """Keeps data in lists so we can capture multiple values eventually, do not know
    what to expect at the moment.

    ChEBI keeps substance names in two places, <NAME> and <SYNONYM> tags.
    """

    etree=ET.parse(filename)
    root=etree.getroot()

    _chebi_id = None
    _definition = []
    _names = set()
    _formula = None
    _smiles = None

    for elem in root:
        if myoptions.debug: print("Level 0: ", elem.tag, ' ', elem.attrib, ' ', elem.text)
        for child in elem:
            if myoptions.debug > 1: print("Items: ", str(child.items()))

            if myoptions.debug: print("Level 1: tag:", child.tag, 'attrib:', child.attrib, 'text:', child.text)
            if child.tag == 'ID':
                if myoptions.debug: print("L1: child.attrib: ", child.attrib, "child.tag: ", child.tag)
                if not _chebi_id:
                    _chebi_id = child.text
                else:
                    raise ValueError("Does ChEBI have multiple IDs for %s ?" % str(child.text))
            if child.tag == 'NAME':
                # ['casbene', 'casbene', 'Casbene']
                # ['terpinolene', 'Terpinolene', 'terpinolene', 'Terpinolen', 'isoterpinene', 'alpha-terpinolene', '4-isopropylidene-1-methylcyclohexene', '1-methyl-4-(1-methylethylidene)cyclohexene', '1-methyl-4-(1-methylethylidene)-1-cyclohexene', '1,4(8)-p-menthadiene']
                # ['epi-cedrol', 'epi-cedrol', '8-epicedrol', '8-epicedrol', '8-epi-cedrol', '(3R,3aS,6S,7R,8aS)-3,6,8,8-tetramethyloctahydro-1H-3a,7-methanoazulen-6-ol', '(-)-epicedrol']
                # ['viridiflorene', 'viridiflorene', 'Ledene', 'Leden']
                # ['dammarenediol-ii', 'dammarenediol-ii', 'dammarenediol ii', 'dammarenediol', 'dammar-24-ene-3beta,20-diol', 'dammar-24-ene-20(ss),3beta-diol', '8-methyl-18-nor-lanost-24-ene-3beta,20-diol', '8-methyl-18-nor-lanost-24-en-3beta,20-diol', '3beta-glucodammar-24-ene-3,20-diol', '(20s)-dammarenediol', '(20s)-dammar-24-ene-3beta,20-diol']
                # ['beta-phellandren', 'beta-phellandrene', '3-isopropyl-6-methylene-1-cyclohexene', '2-p-menthadiene', '3-methylene-6-(1-methylethyl)cyclohexene', '4-isopropyl-1-methylene-2-cyclohexene']
                # ['ophiobolin F', 'ophiobolene']
                # ['(+)-vkiteagnusin d', 'viteagnusin d'] # typo in ChEBI
                _names.update([child.text.lower()])
            if child.tag == 'DEFINITION':
                if not _definition:
                    _definition = [child.text]
                else:
                    _definition += [child.text]
            if child.tag == 'FORMULA':
                if not _formula:
                    _formula = child.text
            if child.tag == 'SYNONYM':
                _names.update([child.text.lower()])
            if child.tag == 'SMILES':
                if not _smiles:
                    _smiles = child.text
            #if child.tag == 'INCHI':
            #    if not _inchi:
            #        _inchi = child.text
    if not _names:
        _names = []
    else:
        _names = list(_names)

    if myoptions.debug: print("Info: IDs: %s, names: %s, definition: %s, formula: %s, smiles: %s" % (str(_chebi_id), str(_names), str(_definition), str(_formula), str(_smiles)))
    return(_chebi_id, list(_names), _definition, _formula, _smiles)


def parse_uniprot_xml(filename, terpenes, uniprot_pri_acc2aliases, uniprot_aliases2pri_acc):
    """Parse a single XML stream (a file pre-fetched into a local cache) from Uniprot.

<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">
  <entry created="2018-11-07" dataset="Swiss-Prot" modified="2019-12-11" version="30">
    <accession>G0LES5</accession>
    <name>TRI5_TRIAR</name>
    ...
  </entry>
  <copyright>
Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms
Distributed under the Creative Commons Attribution (CC BY 4.0) License
</copyright>
</uniprot>
    """

    if not os.path.exists(filename):
        raise ValueError("File %s does not exist." % str(filename))
    else:
        if filename.endswith('xml.gz'):
            _file = gzip.open(filename)
        else:
            _file = open(filename)
    etree = ET.iterparse(_file)

#    try:
#        etree=ET.parse(filename)
#    except ET.ParseError:
#        raise ET.ParseError("Maybe the file %s is not in XML format?" % str(filename))
#    root=etree.getroot() # AttributeError: 'IterParseIterator' object has no attribute 'getroot'

    _accessions = []
    _uniprot_name = None
    _protein_names = [] # sometimes UniProt has unset protein_names and also feature_description
    _feature_descriptions = []
    _chebi_ids = []
    _sequence = None
    _organism = None
    _lineage = None

#    for elem in root:
    for elem_tuple in etree:
        elem = elem_tuple[1]
        # print("Debug: elem=%s" % str(elem))
        if myoptions.debug: print("Level 0: ", elem.tag, ' ', elem.attrib, ' ', elem.text)
        if elem.tag == '{http://uniprot.org/uniprot}entry':
            if _accessions:
                if myoptions.debug > 1: print("Reached items: %s, returning results parsed so far for %s" % (str(elem.items()), _accessions))
                if not _protein_names:
                    # <uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">
                    #   <entry created="2011-10-19" dataset="TrEMBL" modified="2021-04-07" version="68">
                    #     <accession>G1JUH4</accession>
                    #     <name>G1JUH4_SOLLC</name>
                    #     <protein>
                    #       <submittedName>
                    #         <fullName evidence="4 5">Beta myrcene/limonene synthase</fullName>
                    #       </submittedName>
                    raise ValueError("No proteins descriptions were parsed for _accessions=%s" % str(_accessions))
                else:
                    if myoptions.debug: print("Info: Yielding a single entry %s from file %s" % (_accessions[0], str(filename)))
                yield(_accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence)
            _accessions = []
            _uniprot_name = None
            _protein_names = [] # sometimes UniProt has unset protein_names and also feature_description
            _feature_descriptions = []
            _chebi_ids = []
            _sequence = None
            _organism = None
            _lineage = None
        elif elem.tag == '{http://uniprot.org/uniprot}copyright' and _accessions:
            if myoptions.debug > 1: print("Reached items: %s, returning results parsed so far for %s" % (str(elem.items()), _accessions))
            if not _protein_names:
                # <uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">
                #   <entry created="2011-10-19" dataset="TrEMBL" modified="2021-04-07" version="68">
                #     <accession>G1JUH4</accession>
                #     <name>G1JUH4_SOLLC</name>
                #     <protein>
                #       <submittedName>
                #         <fullName evidence="4 5">Beta myrcene/limonene synthase</fullName>
                #       </submittedName>
                raise ValueError("No proteins descriptions were parsed for _accessions=%s" % str(_accessions))
            else:
                if myoptions.debug: print("Info: Yielding the very last entry %s from file %s" % (_accessions[0], str(filename)))
            yield(_accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence)

        for child in elem:
            if myoptions.debug > 1: print("Items: ", str(child.items()))
    
            if myoptions.debug > 1: print("Level 1: tag:", child.tag, 'attrib:', child.attrib, 'text:', child.text)
            if child.tag == '{http://uniprot.org/uniprot}accession':
                if myoptions.debug > 1: print("Info: Came across accession %s" % child.text)
                if myoptions.debug > 1: print("L1: child.attrib: ", child.attrib, "child.tag: ", child.tag)
                if not _accessions:
                    _accessions = [child.text] # "Q6XDB5"
                else:
                    _accessions += [child.text] # "C0PT91"
                    uniprot_aliases2pri_acc[child.text] = _accessions[0] # point to the primary
                    if _accessions[0] not in uniprot_pri_acc2aliases.keys():
                        uniprot_pri_acc2aliases[_accessions[0]] = [child.text]
                    else:
                        uniprot_pri_acc2aliases[_accessions[0]] += [child.text]
            elif child.tag == '{http://uniprot.org/uniprot}name':
                _uniprot_name = child.text # "TPSD2_PICSI"
            elif child.tag == '{http://uniprot.org/uniprot}sequence':
                _sequence = child.text
            if child.tag == '{http://uniprot.org/uniprot}feature':
                if 'type' in child.attrib.keys() and 'description' in child.attrib.keys() and child.attrib['type'] == 'chain':
                    _feature_descriptions += [child.attrib['description']]
                
            for subchild in child:
                if myoptions.debug > 1: print("Level 2: ", subchild.tag, ' ', subchild.attrib, ' ', subchild.text)
                _chebi_ids_local = []
                if child.tag == '{http://uniprot.org/uniprot}organism':
                    if subchild.tag == '{http://uniprot.org/uniprot}name' and subchild.attrib['type'] == 'scientific':
                        _organism = subchild.text
                for sschild in subchild:
                    tag = {}
                    if myoptions.debug > 1: print("Level 3: ", sschild.tag, ' ', sschild.attrib, ' ', sschild.text)
                    if subchild.tag == '{http://uniprot.org/uniprot}recommendedName':
                        if sschild.tag == '{http://uniprot.org/uniprot}fullName':
                            _protein_names = [sschild.text]
                    elif subchild.tag == '{http://uniprot.org/uniprot}alternativeName':
                        if sschild.tag == '{http://uniprot.org/uniprot}fullName':
                            _protein_names += [sschild.text]
                    elif child.tag == '{http://uniprot.org/uniprot}protein':
                        if subchild.tag == '{http://uniprot.org/uniprot}submittedName':
                            if sschild.tag == '{http://uniprot.org/uniprot}fullName':
                                _protein_names += [sschild.text] # G1JUH4
                    elif child.tag == '{http://uniprot.org/uniprot}comment' and 'type' in child.attrib.keys() and child.attrib['type'] == 'catalytic activity' and subchild.tag == '{http://uniprot.org/uniprot}reaction' and sschild.tag == '{http://uniprot.org/uniprot}dbReference' and sschild.attrib['type'] == 'ChEBI':
                        # do not even fetch unwanted ChEBI Ids
                        if sschild.attrib['id'] not in non_terpene_chebi_ids:
                            if sschild.attrib['id'] not in _chebi_ids_local:
                                _chebi_ids_local.append(sschild.attrib['id'])
    
                    # Level 1:  {http://uniprot.org/uniprot}comment   {'type': 'catalytic activity'}   
                    # Level 2:  {http://uniprot.org/uniprot}reaction   {'evidence': '3'}
                    # Level 3:  {http://uniprot.org/uniprot}text   {}   (2E)-geranyl diphosphate = (1S,5S)-alpha-pinene + diphosphate
                    # Level 3:  {http://uniprot.org/uniprot}dbReference   {'id': 'RHEA:25488', 'type': 'Rhea'}   None
                    # Level 3:  {http://uniprot.org/uniprot}dbReference   {'id': 'CHEBI:28660', 'type': 'ChEBI'}   None
                    # Level 3:  {http://uniprot.org/uniprot}dbReference   {'id': 'CHEBI:33019', 'type': 'ChEBI'}   None
                    # Level 3:  {http://uniprot.org/uniprot}dbReference   {'id': 'CHEBI:58057', 'type': 'ChEBI'}   None
                    # Level 3:  {http://uniprot.org/uniprot}dbReference   {'id': '4.2.3.119', 'type': 'EC'}   None
    
                    if child.tag == '{http://uniprot.org/uniprot}organism':
                        if subchild.tag == '{http://uniprot.org/uniprot}lineage':
                            if sschild.tag == '{http://uniprot.org/uniprot}taxon':
                                if not _lineage:
                                    _lineage = [sschild.text]
                                else:
                                    _lineage += [sschild.text]

                if _chebi_ids_local:
                    _chebi_ids += list(set(_chebi_ids_local))
                    _chebi_ids_local = []
    
    if not _protein_names:
        # <uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">
        #   <entry created="2011-10-19" dataset="TrEMBL" modified="2021-04-07" version="68">
        #     <accession>G1JUH4</accession>
        #     <name>G1JUH4_SOLLC</name>
        #     <protein>
        #       <submittedName>
        #         <fullName evidence="4 5">Beta myrcene/limonene synthase</fullName>
        #       </submittedName>
        raise ValueError("No proteins descriptions were parsed for _accessions=%s" % str(_accessions))


def create_cache(cachedir=".TPSdownloader_cache"):
    "make local caching directory for input XML files"

    if not os.path.exists(cachedir):
        os.mkdir(cachedir)
    for subdir in ('uniprot', 'chebi'):
        if not os.path.exists(cachedir + os.path.sep + subdir):
            os.mkdir(cachedir + os.path.sep + subdir)


def downloader(cmdline):
    print("Info:", cmdline)
    _args = shlex_split(cmdline)
    _stdout, _stderr = Popen(_args, shell=False, stdout=PIPE, stderr=PIPE).communicate()
    if _stderr and _stderr[0]:
        sys.stderr.write("Error: '%s' gave '%s'\n" % (cmdline, str(_stderr)))
    sys.stdout.flush()
    sys.stderr.flush()
    time.sleep(30)
    #_handle = curl.Curl(base_url="https://www.uniprot.org/uniprot/" + myid)
    #_myfile = open(path + myid + ".xml", 'wb')
    #_myfile.write(_handle.get())
    #_myfile.close()


def downloader_wrapper(myid, dbname, cachedir, url):
    """curl cannot fetch https://www.uniprot.org/uniprot/Q5Gj59.xml
    but wget can as it accepts redirect to the primary entry. But we
    do not want to store the primary entry under the filename of the
    secondary at least. let's hope we get to the primary entry ID via
    another path.
    """

    if cachedir[-1] != os.path.sep:
        _cachedir = cachedir + os.path.sep + dbname + os.path.sep
    else:
        _cachedir = cachedir + dbname + os.path.sep
    if dbname == 'uniprot' and url[-1] != '/':
        url += '/'
    if os.path.exists(_cachedir):
        _filename = _cachedir + myid + ".xml"
        if os.path.exists(_filename) and not os.path.getsize(_filename):
            os.remove(_filename)
        if not os.path.exists(_filename):
            #print("Debug: fetching %s from uniprot" % myid)
            # older curl version do not support --no-progress-meter
            _cmdline = "curl --silent --show-error -o " + _filename + " " + url + myid
            if dbname == 'uniprot':
                _cmdline += ".xml"
            downloader(_cmdline)
            # curl does not follow redirects if we asked for a secondary/alias accession
            # --2021-05-12 22:49:30--  https://www.uniprot.org/uniprot/D8R8K9.xml
            # Přesměrováno na: /uniprot/G9MAN7.xml [následuji]
            # --2021-05-12 22:49:30--  https://www.uniprot.org/uniprot/G9MAN7.xml
            if os.path.exists(_filename) and not os.path.getsize(_filename):
                sys.stderr.write("Error: Failed to download XML data using '%s' command, re-trying with wget\n" % _cmdline)
                os.remove(_filename)
                _cmdline = "wget --no-proxy --directory-prefix=" + _cachedir + " --max-redirect=1 -o " + _filename + ".log " + url + myid
                if dbname == 'uniprot':
                    _cmdline += ".xml"
                downloader(_cmdline)
            # could also prettyprint the XML files using
            # xml_pp -i.bak _filename
        else:
            if myoptions.debug:
                sys.stderr.write("Debug: File %s already exists\n" % _filename)
    else:
        sys.stderr.write("Error: Directory %s does not exist\n" % _cachedir)


def download_uniprot(myid, path=".TPSdownloader_cache" + os.path.sep + 'uniprot' + os.path.sep):
    """Download a page like
    https://www.uniprot.org/uniprot/A0A2K9RFZ2.xml
    https://www.uniprot.org/uniprot/D8R8K9.xml

    Some entries have multiple Accessions, like

    <accession>Q6XDB5</accession>
    <accession>C0PT91</accession>
    <name>TPSD2_PICSI</name>

    """

    downloader_wrapper(myid, 'uniprot', ".TPSdownloader_cache" + os.path.sep, "https://www.uniprot.org/uniprot/")
    return ".TPSdownloader_cache" + os.path.sep + 'uniprot' + os.path.sep + myid + '.xml'


def download_chebi(myid, path=".TPSdownloader_cache" + os.path.sep + 'chebi' + os.path.sep):
    """Download a page like https://www.ebi.ac.uk/chebi/saveStructure.do?xml=true&chebiId=58622

    ChEBI also provides SQL dumps, probably an overkill for our
    purpose.
    """

    downloader_wrapper(myid, 'chebi', ".TPSdownloader_cache" + os.path.sep, "https://www.ebi.ac.uk/chebi/saveStructure.do?xml=true&chebiId=")


def substance_of_interest():
    """Evaluate structure of chemical substances involved in a chemcical
    reaction. If they are not water, GGPP, FFPP, return True so we can output
    them.

    Try just checking if the substance contains more than 5 carbon atoms.
    Currently we get around this by blacklisting some ChEBI Id's we do not need.
    """

    pass


def print_df(df):
    """
    label=='Uniprot ID'
    label=='Name'
    """

    for label, content in df.items():
        print(f'label: {label}')
        print(f'content: {content}', sep='\n')


def write_csv_and_xls(df):
    """Write CSV output, per one Uniprot ID A0A2K9RFZ2 output even
    multiple lines if there are multiple reactions catalyzed

    https://www.uniprot.org/uniprot/A0A2K9RFZ2
    https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:58622
    """
    
    _datetime = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    df.to_csv("TPSdownloader_" + _datetime + ".csv", index=False)
    df.to_excel("TPSdownloader_" + _datetime + ".xls", sheet_name="Sheet1", index=False)


def sanitize_input_text_values(i):
    "Convert string representation of pythonic values into real python datatypes."

    if i is np.NaN:
        _i = None
    elif i is None:
        _i = None
    else:
        _i = i

    if i and isinstance(i, bool):
        _i = True
    elif not i and isinstance(i, bool):
        _i = None
    elif i and isinstance(i, str):
        _i = i.strip()

    if _i is np.NaN:
        _i = None
    elif _i is None:
        pass
    elif _i == 'nan':
        _i = None
    elif str(_i) == 'nan':
        _i = None

    if not _i:
        return None

    if _i.startswith('[') and _i.endswith(']'):
        _sub_i = parse_list_or_set(_i) # split the list text representation into entries
    else:
        _sub_i = [_i]
    return _sub_i


def parse_storage(filename):
    """Parse XLS into memory, into pythonic lists and sets.
    """

    _uniprot_dict_of_lists, _chebi_dict_of_lists, _copy_without_chebi_id, _empty_template_dict_of_lists = initialize_data_structures()
    _storage_df = pd.read_excel(filename, 'Sheet1', index_col=None, na_values=["NA"])
    for _colname in ['Uniprot ID', 'ChEBI ID', 'Name', 'Description', 'Species', 'Taxonomy', 'Amino acid sequence', 'Kingdom (plant, fungi, bacteria)', 'Notes', 'Publication (URL)', 'Substrate ChEBI IDs', 'Intermediate ChEBI IDs', 'Product ChEBI IDs', 'Cofactors ChEBI IDs'] + extra_product_colnames + extra_substrate_colnames + extra_cofactor_colnames + extra_intermediate_colnames:
        for _val in _storage_df[_colname]:
            _uniprot_dict_of_lists[_colname].append(sanitize_input_text_values(_val))

    for _colname in ['ChEBI ID', 'Compound name', 'Compound description', 'Formula', 'SMILES', 'Type (mono, sesq, di, …)']:
        for _val in _storage_df[_colname]:
            if _val not in _chebi_dict_of_lists[_colname].values(): # remove redundancies
                _chebi_dict_of_lists[_colname].append(sanitize_input_text_values(_val))

    del(_storage_df)
    return _uniprot_dict_of_lists, _chebi_dict_of_lists, _copy_without_chebi_id, _empty_template_dict_of_lists


def process_parsed_uniprot_values(all_uniprot_ids, all_chebi_ids, uniprot_dict_of_lists, chebi_dict_of_lists, already_parsed, accessions, chebi_ids, protein_names, feature_descriptions, organism, lineage, sequence, uniprot_pri_acc2aliases, uniprot_aliases2pri_acc):
    if accessions[0] not in already_parsed:
        print("Info: Storing uniprot values for entry %s into uniprot_dict_of_lists" % accessions[0])
        for _uniprot_id in accessions:
            all_uniprot_ids.update([_uniprot_id])
        for _chebi_id in chebi_ids:
            all_chebi_ids.update([_chebi_id])
        if len(accessions) > 1:
            if accessions[0] not in uniprot_pri_acc2aliases.values():
                for _alias in accessions[1:]:
                    uniprot_aliases2pri_acc[_alias] = accessions[0]
                uniprot_pri_acc2aliases[accessions[0]] = accessions[1:]
            else:
                if uniprot_aliases2pri_acc[_alias] != accessions[0]:
                    raise ValueError("Primary accession %s already recorded is not same as the parsed right now: %s" % (uniprot_aliases2pri_acc[_alias], accessions[0]))

        uniprot_dict_of_lists['Uniprot ID'].append(accessions)
        uniprot_dict_of_lists['ChEBI ID'].append(chebi_ids)

        # parse values for ChEBI entries mentioned in the UniProt record
        _substrate_ids = []
        _intermediate_ids = []
        _product_ids = []
        _cofactor_ids = []
        _old_vector = (False, False, False, False)
        _has_product, _has_substrate, _has_intermediate, _has_cofactor = None, None, None, None
        _terpene_type_saved = ''
        for _chebi_id in chebi_ids:
            download_chebi(_chebi_id)
            _filename = ".TPSdownloader_cache" + os.path.sep + 'chebi' + os.path.sep + _chebi_id + '.xml'
            if os.path.exists(_filename) and os.path.getsize(_filename):
                _chebi_id2, _names, _definition, _formula, _smiles = parse_chebi_xml(_filename)
                if _formula:
                    _terpene_type = classify_terpene(_formula)
                else:
                    _terpene_type = None
    
                _has_product, _has_substrate, _has_intermediate, _has_cofactor = split_chebi_data_into_substrates_and_products(_chebi_id2, _terpene_type, _substrate_ids, _intermediate_ids, _product_ids, _cofactor_ids)
                _old_vector = [x | y for x, y in zip(tuple(_old_vector), (_has_product, _has_substrate, _has_intermediate, _has_cofactor))] # perform bitwise OR operation to merge the bits
            else:
                _has_product, _has_substrate, _has_intermediate, _has_cofactor = None, None, None, None
                _terpene_type = None

            if _terpene_type:
                _terpene_type_saved = str(_terpene_type)

            if _chebi_id2 and _chebi_id2 not in chebi_dict_of_lists['ChEBI ID']:
                chebi_dict_of_lists['ChEBI ID'].append(_chebi_id2)
                chebi_dict_of_lists['Compound name'].append(_names)
                chebi_dict_of_lists['Compound description'].append(_definition)
                chebi_dict_of_lists['Formula'].append(_formula)
                chebi_dict_of_lists['SMILES'].append(_smiles)
                if _terpene_type:
                    chebi_dict_of_lists['Type (mono, sesq, di, …)'].append(_terpene_type)
                else:
                    chebi_dict_of_lists['Type (mono, sesq, di, …)'].append('')

        if _has_substrate:
            uniprot_dict_of_lists['Substrate ChEBI IDs'].append(_substrate_ids)
        else:
            uniprot_dict_of_lists['Substrate ChEBI IDs'].append([])
        if _has_intermediate:
            uniprot_dict_of_lists['Intermediate ChEBI IDs'].append(_intermediate_ids)
        else:
            uniprot_dict_of_lists['Intermediate ChEBI IDs'].append([])
        if _has_product:
            uniprot_dict_of_lists['Product ChEBI IDs'].append(_product_ids)
        else:
            uniprot_dict_of_lists['Product ChEBI IDs'].append([])
        if _has_cofactor:
            uniprot_dict_of_lists['Cofactors ChEBI IDs'].append(_cofactor_ids)
        else:
            uniprot_dict_of_lists['Cofactors ChEBI IDs'].append([])

        #uniprot_dict_of_lists['ChEBI ID'].append(_substrate_ids + _intermediate_ids + _product_ids + _cofactor_ids)

        if not _has_product:
            for _colname in extra_product_colnames:
                uniprot_dict_of_lists[_colname].append([])
        else:
            _indexes = [chebi_dict_of_lists['ChEBI ID'].index(x) for x in _product_ids]
            #_uniprot_index = uniprot_dict_of_lists['ChEBI ID'].index(accessions[0])
            uniprot_dict_of_lists['Name of product'].append([chebi_dict_of_lists['Compound name'][x] for x in _indexes])
            uniprot_dict_of_lists['Product compound description'].append([chebi_dict_of_lists['Compound description'][x] for x in _indexes])
            uniprot_dict_of_lists['Chemical formula of product'].append([chebi_dict_of_lists['Formula'][x] for x in _indexes])
            uniprot_dict_of_lists['SMILES of product (including stereochemistry)'].append([chebi_dict_of_lists['SMILES'][x] for x in _indexes])

        if not _has_substrate:
            for _colname in extra_substrate_colnames:
                uniprot_dict_of_lists[_colname].append([])
        else:
            _indexes = [chebi_dict_of_lists['ChEBI ID'].index(x) for x in _substrate_ids]
            uniprot_dict_of_lists['Substrate (including stereochemistry)'].append([chebi_dict_of_lists['Compound name'][x] for x in _indexes])
            uniprot_dict_of_lists['Substrate compound description'].append([chebi_dict_of_lists['Compound description'][x] for x in _indexes])
            uniprot_dict_of_lists['Chemical formula of substrate'].append([chebi_dict_of_lists['Formula'][x] for x in _indexes])
            uniprot_dict_of_lists['SMILES of substrate (including stereochemistry)'].append([chebi_dict_of_lists['SMILES'][x] for x in _indexes])

        if not _has_intermediate:
            for _colname in extra_intermediate_colnames:
                uniprot_dict_of_lists[_colname].append([])
        else:
            _indexes = [chebi_dict_of_lists['ChEBI ID'].index(x) for x in _intermediate_ids]
            uniprot_dict_of_lists['Name of intermediate'].append([chebi_dict_of_lists['Compound name'][x] for x in _indexes])
            uniprot_dict_of_lists['Intermediate compound description'].append([chebi_dict_of_lists['Compound description'][x] for x in _indexes])
            uniprot_dict_of_lists['Chemical formula of intermediate'].append([chebi_dict_of_lists['Formula'][x] for x in _indexes])
            uniprot_dict_of_lists['SMILES of intermediate (including stereochemistry)'].append([chebi_dict_of_lists['SMILES'][x] for x in _indexes])

        if not _has_cofactor:
            for _colname in extra_cofactor_colnames:
                uniprot_dict_of_lists[_colname].append([])
        else:
            _indexes = [chebi_dict_of_lists['ChEBI ID'].index(x) for x in _cofactor_ids]
            uniprot_dict_of_lists['Cofactors'].append([chebi_dict_of_lists['Compound name'][x] or x in _indexes])
            uniprot_dict_of_lists['Cofactors compound description'].append([chebi_dict_of_lists['Compound description'][x] for x in _indexes])
            uniprot_dict_of_lists['Chemical formula of cofactor'].append([chebi_dict_of_lists['Formula'][x] for x in _indexes])
            uniprot_dict_of_lists['SMILES of cofactor (including stereochemistry)'].append([chebi_dict_of_lists['SMILES'][x] for x in _indexes])

        uniprot_dict_of_lists['Name'].append(protein_names)
        uniprot_dict_of_lists['Description'].append(feature_descriptions)
        uniprot_dict_of_lists['Species'].append(organism)
        uniprot_dict_of_lists['Taxonomy'].append(lineage)
        uniprot_dict_of_lists['Amino acid sequence'].append(sequence)
        if not lineage:
            _kingdom = ''
        elif 'Bacteria' == lineage[0]: # some taxons are only assigned as ['Bacteria'], without further specs
            _kingdom = 'Bacteria'
        elif 'Viridiplantae' == lineage[1]:
            _kingdom = 'Plantae'
        elif 'Fungi' == lineage[1]:
            _kingdom = 'Fungi'
        elif 'Homo sapiens' in lineage:
            _kingdom = 'Human'
        elif 'Animalia' in lineage:
            _kingdom = 'Animal'
        else:
            _kingdom = 'unknown'
        uniprot_dict_of_lists['Kingdom (plant, fungi, bacteria)'].append(_kingdom)
        uniprot_dict_of_lists['Notes'].append('')
        uniprot_dict_of_lists['Publication (URL)'].append('')
        already_parsed.append(accessions[0])
    else:
        print("Info: Accession %s already kept in uniprot_dict_of_lists" % accessions[0])

    for x in uniprot_dict_of_lists.keys():
        print("Info: uniprot_dict_of_lists[%s]=%d" % (x, len(uniprot_dict_of_lists[x])))


def fetch_ids_from_xlsx(filename, terpenes, uniprot_pri_acc2aliases, uniprot_aliases2pri_acc, uniprot_dict_of_lists, already_parsed, all_uniprot_ids, all_chebi_ids):
    _ids = []
    _df = pd.read_excel(filename, 'Sheet1', index_col=None, na_values=["NA"])
    for i in _df['Uniprot ID']:
        _sub_i = sanitize_input_text_values(i)

        if len(_sub_i) > 1:
            for _ii in _sub_i:
                print("Looping over %s from %s, originates from %s" % (_ii, filename, str(i)))
                _id = None
                # remove redundancies but keep ordering
                if _ii is not None and _ii not in already_parsed:
                    # check if this is a primary accession, if not, convert it to primary
                    if _ii in uniprot_pri_acc2aliases.keys():
                        # already known primary accession, probably already in cache but maybe inferred from other sources
                        _filename = download_uniprot(_ii)
                        _id = _ii
                    elif _ii in uniprot_pri_acc2aliases.values():
                        # is a secondary/alias accession
                        _iii = uniprot_aliases2pri_acc[_ii]
                        _id = _iii
                        if _iii in already_parsed:
                            # we already have parsed the primary acc uniprot entry into memory
                            # for safety re-try fetching the data if it is not in the cache
                            print("Info: Entry %s is an alias pointing to %s which we already have, skipping download" % (_ii, _iii))
                            # _filename = download_uniprot(_id)
                            _filename = None # prevent duplicated parsing of the input XML
                        elif _iii not in uniprot_aliases2pri_acc.keys():
                            print("Info: Entry %s is an alias pointing to %s which we already have, also skipping download" % (_ii, _iii))
                            _filename = download_uniprot(_iii)
                            #_filename = None
                    else:
                        # new or already known primary accession
                        print("Info: Entry %s is new primary accession, downloading if not existing yet" % (_ii))
                        _filename = download_uniprot(_ii)
                        _id = _ii
                if _id:
                    _ids.append(_id)
        else:
            # there are no secondary accessions
            _ids += _sub_i
            if _sub_i[0] not in already_parsed:
                _filename = download_uniprot(_sub_i[0])
    return set(_ids)


def parse_multi_entry_uniprot_xml(filename, terpenes, uniprot_pri_acc2aliases, uniprot_aliases2pri_acc, uniprot_dict_of_lists, already_parsed):
    if not os.path.exists(filename):
        raise ValueError("File %s does not exist." % str(filename))
    else:
        if filename.endswith('xml.gz'):
            _file = gzip.open(filename)
        else:
            _file = open(filename)
    _xml_tree = ET.iterparse(_file)
    print("Info: ET.iterparse(%s) returned %s" % (str(_file), str(_xml_tree)))
    for _items in _xml_tree:
        print("Info: Acting upon %s" % str(_items))
        for _uniprot_entry in _items:
            print("Info: Acting over entry %s in %s" % (str(_uniprot_entry), filename))
            # <Element '{http://uniprot.org/uniprot}entry'
            (_accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence) = parse_uniprot_xml(_uniprot_entry, terpenes, uniprot_pri_acc2aliases, uniprot_aliases2pri_acc)
            process_parsed_uniprot_values(_uniprot_id, _chebi_id, _all_uniprot_ids, _all_chebi_ids, dict_of_lists, already_parsed, _accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence, uniprot_pri_acc2aliases, uniprot_aliases2pri_acc)


_r1 = re.compile(r'C[0-9]+')

def classify_terpene(formula):
    _match = _r1.search(formula)
    if _match:
        _carbon_count = int(formula[_match.span()[0]+1:_match.span()[1]])
    else:
        # CHEBI:35757
        # <FORMULA>CO2R</FORMULA>
        _carbon_count = 0 # is zero or one but we do not care in this case anyway
    if _carbon_count > 12 and _carbon_count < 17:
        _terpene_type = 'sesq'
    elif _carbon_count > 8 and _carbon_count < 13:
        _terpene_type = 'mono' # C12
    elif _carbon_count > 18 and _carbon_count < 22:
        _terpene_type = 'di'
    elif _carbon_count > 28 and _carbon_count < 32:
        _terpene_type = 'tri'
    elif _carbon_count > 31 and _carbon_count < 37:
        _terpene_type = 'sesquar'
    elif _carbon_count > 38 and _carbon_count < 42:
        _terpene_type = 'tetra'
    elif _carbon_count > 23 and _carbon_count < 27:
        _terpene_type = 'sest'
    elif not _carbon_count:
        # CHEBI:35757
        # <FORMULA>CO2R</FORMULA>
        _terpene_type = None
    else:
        _terpene_type = 'unexpected'
    return _terpene_type


def initialize_data_structures():
    _uniprot_dict_of_lists = {'Uniprot ID': [], 'Name': [], 'Description': [], 'Species': [], 'Taxonomy': [], 'Amino acid sequence': [], 'Kingdom (plant, fungi, bacteria)': [], 'ChEBI ID': [], 'Substrate ChEBI IDs': [], 'Intermediate ChEBI IDs': [], 'Product ChEBI IDs': [], 'Cofactors ChEBI IDs': [], 'Notes': [], 'Publication (URL)': []}

    _chebi_dict_of_lists = {'ChEBI ID': [], 'Compound name': [], 'Compound description': [], 'Formula': [], 'SMILES': [], 'Type (mono, sesq, di, …)': []}

    for _colname in extra_product_colnames + extra_substrate_colnames + extra_cofactor_colnames + extra_intermediate_colnames:
        _uniprot_dict_of_lists[_colname] = []

    _copy_without_chebi_id = copy.deepcopy(_chebi_dict_of_lists)
    _copy_without_chebi_id.pop('ChEBI ID')
    _empty_template_dict_of_lists = copy.deepcopy(_uniprot_dict_of_lists)
    _empty_template_dict_of_lists.update(_copy_without_chebi_id)
    return _uniprot_dict_of_lists, _chebi_dict_of_lists, _copy_without_chebi_id, _empty_template_dict_of_lists


def split_chebi_data_into_substrates_and_products(_chebi_id2, _terpene_type, substrate_ids, intermediate_ids, product_ids, cofactor_ids):
    "This function needs reworking to simultaneously split into all groups."

    _has_product = False
    _has_substrate = False
    _has_intermediate = False
    _has_cofactor = False

    if _terpene_type:
        product_ids.append(_chebi_id2)
        _has_product = True
    elif _chebi_id2 in ['CHEBI:17211', 'CHEBI:14299', 'CHEBI:5332', 'CHEBI:42877', 'CHEBI:24223', 'CHEBI:58635', 'CHEBI:30939', 'CHEBI:10760', 'CHEBI:29558', 'CHEBI:57907', 'CHEBI:58756', 'CHEBI:58057', 'CHEBI:162247', 'CHEBI:60374', 'CHEBI:138307', 'CHEBI:61984', 'CHEBI:64801', 'CHEBI:15441', 'CHEBI:18728', 'CHEBI:11026', 'CHEBI:11072', 'CHEBI:372', 'CHEBI:58622', 'CHEBI:58206', 'CHEBI:64283', 'CHEBI:138223', 'CHEBI:58553', 'CHEBI:57533', 'CHEBI:7525', 'CHEBI:138305', 'CHEBI:15440', 'CHEBI:10843', 'CHEBI:9245', 'CHEBI:10795', 'CHEBI:15104', 'CHEBI:26746', 'CHEBI:64283', 'CHEBI:175763', 'CHEBI:17407', 'CHEBI:12874', 'CHEBI:11491', 'CHEBI:42496', 'CHEBI:10700', 'CHEBI:11488', 'CHEBI:12854', 'CHEBI:19789', 'CHEBI:138890', 'CHEBI:138232']: # substrates
        # CHEBI:17211 CHEBI:14299, CHEBI:5332, CHEBI:42877, CHEBI:24223 GPP
        # CHEBI:58635 nebo CHEBI:30939 CHEBI:10760, CHEBI:29558 (+)-copalyl diphosphate
        # CHEBI:57907 (2E,6E,10E,14E)-GFPP
        # CHEBI:58756 (2E,6E,10E)-GGPP
        # CHEBI:58057 (2E)-GPP
        # CHEBI:162247 (2Z,6E)-FPP
        # CHEBI:60374 (2Z,6Z)-FPP
        # CHEBI:138307 (3S,22S)-2,3:22,23-diepoxy-2,3,22,23-tetrahydrosqualene
        # CHEBI:61984 (E)-2-MeGPP alias (E)-2-methylgeranyl diphosphate
        # CHEBI:64801 (R)-tetraprenyl-β-curcumene
        # CHEBI:15441 CHEBI:18728, CHEBI:11026, CHEBI:11072, CHEBI:372 (S)-2,3-epoxysqualene
        # CHEBI:58622 9α-copalyl PP
        # CHEBI:58206 all-trans-heptaprenyl PP
        # CHEBI:64283 copal-8-ol diphosphate(3−)
        # CHEBI:138223 ent-copal-8-ol diphosphate(3−)
        # CHEBI:58553 ent-copalyl diphosphate
        # CHEBI:57533 GGPP
        # CHEBI:7525 NPP
        # CHEBI:138305 pre-α-onocerin
        # CHEBI:15440 CHEBI:10843, CHEBI:9245, CHEBI:10795, CHEBI:15104, CHEBI:26746 squalene
        # CHEBI:64283 8-hydroxycopalyl diphosphate
        # (2E,6E)-FPP(3-) CHEBI:175763
        # CHEBI:17407 (E,E)-FPP FPP CHEBI:12874, CHEBI:11491, CHEBI:42496, CHEBI:10700, CHEBI:11488, CHEBI:12854, CHEBI:19789 (2-trans,6-trans-farnesyl diphosphate) alias (E,E)-FPP alias (2E,6E)-FPP
        # peregrinol diphosphate CHEBI:138890
        # peregrinol PP alias? peregrinol diphosphate(3−) CHEBI:138232
        # 
        # Nevyreseno:
        # (Z,Z)-FPP
        # NNPP
        # 
        substrate_ids.append(_chebi_id2)
        _has_substrate = True
    elif _chebi_id2 in ['CHEBI:18420', 'CHEBI:15377', 'CHEBI:29035']: # cofactors
        # CHEBI:18420 Mg2+
        # CHEBI:29035 Mn2+
        # CHEBI:15377 H2O
        cofactor_ids.append(_chebi_id2)
        _has_cofactor = True
    elif ['CHEBI:63190', 'CHEBI:58622', 'CHEBI:63190', 'CHEBI:58553', 'CHEBI:64283', 'CHEBI:58635', 'CHEBI:30939', 'CHEBI:10760', 'CHEBI:29558']:
        # CHEBI:63190 (+)-β-caryophyllene
        # CHEBI:58622 9α-copalyl diphosphate
        # CHEBI:63190 (S)-β-bisabolene
        # CHEBI:58553 ent-copalyl diphosphate
        # CHEBI:64283 copal-8-ol diphosphate(3−)
        # CHEBI:58635 CHEBI:30939 CHEBI:10760, CHEBI:29558 (+)-copalyl diphosphate
        intermediate_ids.append(_chebi_id2)
        _has_intermediate = True
    else:
        raise ValueError("Unexpected compound %s" % _chebi_id2)

    return _has_product, _has_substrate, _has_intermediate, _has_cofactor

        

def print_dict_lengths(somedict, dictname):
    for _item in somedict.keys():
        if myoptions.debug > 3: print("%s: Key=%s, len=%d" % (dictname, _item, len(somedict[_item])))
        print("Info: Will output into CSV and XLS files in total %s: Key=%s, len=%d" % (dictname, _item, len(somedict[_item])))

def main():
    create_cache()
    _terpenes = parse_known_terpenes()
    _uniprot_pri_acc2aliases = {}
    _uniprot_aliases2pri_acc = {}

    _all_uniprot_ids = set()
    _all_chebi_ids = set()

    if myoptions.xls_storage and os.path.exists(myoptions.xls_storage) and myoptions.xls_storage != 'None':
        _uniprot_dict_of_lists, _chebi_dict_of_lists, _copy_without_chebi_id, _empty_template_dict_of_lists = parse_storage(myoptions.xls_storage)
    else:
        # one can disable the default value with passing None on the commandline
        _uniprot_dict_of_lists, _chebi_dict_of_lists, _copy_without_chebi_id, _empty_template_dict_of_lists = initialize_data_structures()

    _output_dict_of_lists = copy.deepcopy(_empty_template_dict_of_lists)

    if myoptions.debug: print("Debug: Initialized data structures and parsed XLS storage into _uniprot_dict_of_lists=%s", str(_uniprot_dict_of_lists))

    _already_parsed = []
    #  parse previously obtained multi-entry XML data, if any
    for _filename in os.listdir('.TPSdownloader_cache/uniprot/multientry/'):
        print("Info: Found multi-entry XML file %s" % '.TPSdownloader_cache/uniprot/multientry/' + _filename)
        if os.path.getsize('.TPSdownloader_cache/uniprot/multientry/' + _filename):
            print("Info: Parsing %s" % '.TPSdownloader_cache/uniprot/multientry/' + _filename)
            for _accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence in parse_uniprot_xml('.TPSdownloader_cache/uniprot/multientry/' + _filename, _terpenes, _uniprot_pri_acc2aliases, _uniprot_aliases2pri_acc):
                process_parsed_uniprot_values(_all_uniprot_ids, _all_chebi_ids, _uniprot_dict_of_lists, _chebi_dict_of_lists, _already_parsed, _accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence, _uniprot_pri_acc2aliases, _uniprot_aliases2pri_acc)

    if myoptions.debug: print("Info: _all_uniprot_ids=%s" % str(_all_uniprot_ids))
    if myoptions.debug: print("Info: _all_chebi_ids=%s" % str(_all_chebi_ids))

    if myoptions.uniprot_id:
        download_uniprot(myoptions.uniprot_id)
    elif myoptions.uniprot_ids_from_file and os.path.exists(myoptions.uniprot_ids_from_file):
        # get list of accessions, fetch their single-entry XML files unless already in local cache and parse them
        _ids = fetch_ids_from_xlsx(myoptions.uniprot_ids_from_file, _terpenes, _uniprot_pri_acc2aliases, _uniprot_aliases2pri_acc, _uniprot_dict_of_lists, _already_parsed, _all_uniprot_ids, _all_chebi_ids)

    _aliases = _uniprot_aliases2pri_acc.keys()
    for _id in _ids:
        _filename = '.TPSdownloader_cache/uniprot/' + _id + '.xml'
        if _id not in _already_parsed and _id not in _aliases:
            if os.path.exists(_filename) and os.path.getsize(_filename):
                for _accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence in parse_uniprot_xml(_filename, _terpenes, _uniprot_pri_acc2aliases, _uniprot_aliases2pri_acc):
                    process_parsed_uniprot_values(_all_uniprot_ids, _all_chebi_ids, _uniprot_dict_of_lists, _chebi_dict_of_lists, _already_parsed, _accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence, _uniprot_pri_acc2aliases, _uniprot_aliases2pri_acc)
            else:
                print("Info: No such file %s" % str(_filename))
        else:
            print("Info: Skipping %s which is already in _all_uniprot_ids, supposedly a secondary accession" % _id)

    if myoptions.debug: print("Debug: _all_uniprot_ids=%s" % str(_all_uniprot_ids))
    if myoptions.debug: print("Debug: _all_chebi_ids=%s" % str(_all_chebi_ids))

    print("Info: len(_already_parsed)=%s" % len(_already_parsed))
    #print("AAAAA: len(_uniprot_dict_of_lists)=%d, len(_uniprot_dict_of_lists['Uniprot ID'])=%s, _uniprot_dict_of_lists: %s" % (len(_uniprot_dict_of_lists), len(_uniprot_dict_of_lists['Uniprot ID']), str(_uniprot_dict_of_lists)))

    if myoptions.debug > 1: print("Debug: After parsing ChEBI files _uniprot_dict_of_lists=", str(_uniprot_dict_of_lists))
    print("Info: There are %s 'ChEBI ID' entries in _chebi_dict_of_lists: %s" % (len(_chebi_dict_of_lists['ChEBI ID']), str(_chebi_dict_of_lists)))

    # re-copy the parsed data into rows if there are multiple _chebi_ids annotated (pointing to multiple cyclic terpenes), other ChEBI IDs were discarded
    _unique_uniprot_ids = list(_uniprot_dict_of_lists['Uniprot ID'])
    for _uniprot_row_pos, _uniprot_entry in enumerate(_unique_uniprot_ids):
        _chebi_ids = _uniprot_dict_of_lists['ChEBI ID'][_uniprot_row_pos]
        # _chebi_ids = map(lambda x: x[1:], _uniprot_dict_of_lists['ChEBI ID'][_uniprot_row_pos])
        # for each CHEBI item in _chebi_ids, output a dedicated line in the output
        if len(_output_dict_of_lists['ChEBI ID']) != len(_output_dict_of_lists['Uniprot ID']):
            print_dict_lengths(_uniprot_dict_of_lists, '_uniprot_dict_of_lists')
            print_dict_lengths(_chebi_dict_of_lists, '_chebi_dict_of_lists')
            print_dict_lengths(_output_dict_of_lists, '_output_dict_of_lists')
            raise ValueError("Sizes do not match,\n_output_dict_of_lists: %s\n" % str(_output_dict_of_lists))
        if _chebi_ids and list(_chebi_ids):
            for _chebi_ids_index, _chebi_id in enumerate(_chebi_ids):
                _chebi_row_pos = _chebi_dict_of_lists['ChEBI ID'].index(_chebi_id)
                if _chebi_id:
                    _chebi_id2 = _chebi_dict_of_lists['ChEBI ID'][_chebi_row_pos]
                else:
                    _chebi_id2 = None

                if _chebi_id and _chebi_id != set() and _chebi_id != _chebi_id2:
                    sys.stderr.write("Error: _chebi_id=%s != _chebi_id2=%s\n" % (str(_chebi_id), str(_chebi_id2)))
                if myoptions.debug > 1: print("Found ChEBI ID in Uniprot data: %s at row %s, also in ChEBI data: %s at row %s" % (_chebi_id, _uniprot_row_pos, _chebi_id2, _chebi_row_pos))
                if myoptions.debug > 1: print("  _chebi_row_pos=%s" % _chebi_row_pos)

                if _chebi_id == _chebi_ids[_chebi_ids_index]:
                    # re-copy the Uniprot-originating data
                    for _column in _uniprot_dict_of_lists.keys():
                        if _column != 'ChEBI ID':
                            if _column in ['Substrate ChEBI IDs', 'Intermediate ChEBI IDs', 'Product ChEBI IDs', 'Cofactors ChEBI IDs'] + extra_product_colnames + extra_substrate_colnames + extra_cofactor_colnames + extra_intermediate_colnames:
                                try:
                                    _val = _uniprot_dict_of_lists[_column][_uniprot_row_pos][_chebi_ids_index] # slice out for 'Product ChEBI IDs', 'Name of product', 'Product compound description', 'Chemical formula of product', 'SMILES of product (including stereochemistry)' columns using _chebi_ids_index
                                except IndexError:
                                    _val = _uniprot_dict_of_lists[_column][_uniprot_row_pos] # well, some column are still totally empty
                            else:
                                _val = _uniprot_dict_of_lists[_column][_uniprot_row_pos]
                            if _val:
                                _output_dict_of_lists[_column].append(_val)
                            else:
                                _output_dict_of_lists[_column].append('')

                    # fill-in the ChEBI-originating data
                    # copy-out values from _chebi_dict_of_lists but only if they match ChEBI ID stored in _uniprot_dict_of_lists['Uniprot ID'] or are on the same line in _chebi_dict_of_lists
                    # if _uniprot_entry[0] == 'A0A1D6LTV0':
                    #     _debug = 1
                    # else:
                    #     _debug = 0

                    # Uniprot ID: A0A1D6LTV0
                    # _chebi_dict_of_lists[_column][_chebi_row_pos] = CHEBI:128
                    # _uniprot_dict_of_lists['ChEBI ID'][_uniprot_row_pos] = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884']
                    # _chebi_ids = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884'], _chebi_id2 = CHEBI:128, _indexes = [0]
                    # _chebi_dict_of_lists[_column][_chebi_row_pos] = CHEBI:15383
                    # _uniprot_dict_of_lists['ChEBI ID'][_uniprot_row_pos] = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884']
                    # _chebi_ids = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884'], _chebi_id2 = CHEBI:15383, _indexes = [0]
                    # _chebi_dict_of_lists[_column][_chebi_row_pos] = CHEBI:10577
                    # _uniprot_dict_of_lists['ChEBI ID'][_uniprot_row_pos] = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884']
                    # _chebi_ids = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884'], _chebi_id2 = CHEBI:10577, _indexes = [0]
                    # _chebi_dict_of_lists[_column][_chebi_row_pos] = CHEBI:9457
                    # _uniprot_dict_of_lists['ChEBI ID'][_uniprot_row_pos] = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884']
                    # _chebi_ids = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884'], _chebi_id2 = CHEBI:9457, _indexes = [0]
                    # _chebi_dict_of_lists[_column][_chebi_row_pos] = CHEBI:78884
                    # _uniprot_dict_of_lists['ChEBI ID'][_uniprot_row_pos] = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884']
                    # _chebi_ids = ['CHEBI:128', 'CHEBI:15383', 'CHEBI:10577', 'CHEBI:9457', 'CHEBI:78884'], _chebi_id2 = CHEBI:78884, _indexes = [0]
                    for _column in _chebi_dict_of_lists.keys():
                        if _column in ['ChEBI ID', 'Substrate ChEBI IDs', 'Intermediate ChEBI IDs', 'Product ChEBI IDs', 'Cofactors ChEBI IDs']:
                            #if _debug:
                            #    print("_chebi_dict_of_lists[_column][_chebi_row_pos] = %s" % str(_chebi_dict_of_lists[_column][_chebi_row_pos]))
                            #    print("_uniprot_dict_of_lists['ChEBI ID'][_uniprot_row_pos] = %s" % str(_uniprot_dict_of_lists['ChEBI ID'][_uniprot_row_pos]))
                            #    print("_chebi_ids = %s, _chebi_id2 = %s" % (str(_chebi_ids), str(_chebi_id2)))
                            #    print("_chebi_ids_index=%s, _chebi_row_pos=%s" % (_chebi_ids_index, _chebi_row_pos))
                            #    print("Column names: %s" % str(_chebi_dict_of_lists.keys()))
                            _val = _chebi_dict_of_lists[_column][_chebi_row_pos]
                            if _val:
                                _output_dict_of_lists[_column].append(_val)
                            else:
                                _output_dict_of_lists[_column].append('')
                        else:
                            # _chebi_dict_of_lists column names: ['ChEBI ID', 'Compound name', 'Compound description', 'Formula', 'SMILES', 'Type (mono, sesq, di, …)']
                            # _output_dict_of_lists column names: ['Product ChEBI IDs', 'Name of product', 'Product compound description', 'Chemical formula of product', 'SMILES of product (including stereochemistry)']
                            if _column in ['Compound name', 'Compound description', 'Formula', 'SMILES']:
                                try:
                                    _val = _chebi_dict_of_lists[_column][_chebi_row_pos][_chebi_ids_index] # slice out for 'Product ChEBI IDs', 'Name of product', 'Product compound description', 'Chemical formula of product', 'SMILES of product (including stereochemistry)' columns using _chebi_ids_index
                                except IndexError:
                                    _val = _chebi_dict_of_lists[_column][_chebi_row_pos]
                                    # raise IndexError("Cannot slice over column %s" % _column) # well, some column are still totally empty
                            else:
                                _val = _chebi_dict_of_lists[_column][_chebi_row_pos]
                            if _val:
                                _output_dict_of_lists[_column].append(_val)
                            else:
                                _output_dict_of_lists[_column].append('')
        else:
            # re-copy just the Uniprot-originating data
            for _column in _uniprot_dict_of_lists.keys():
                if _column != 'ChEBI ID':
                    try:
                        _val = _uniprot_dict_of_lists[_column][_uniprot_row_pos]
                    except IndexError:
                        sys.stderr.write("Error: There are the following columns defined in _uniprot_dict_of_lists: %s\n" % str(_uniprot_dict_of_lists.keys()))
                        raise IndexError("Row %s in '%s' column is missing, cannot copy its values from _uniprot_dict_of_lists" % (str(_uniprot_row_pos), str(_column)))
                    if _val:
                        _output_dict_of_lists[_column].append(_val)
                    else:
                        _output_dict_of_lists[_column].append('')

            # fill-in the ChEBI-originating data
            for _column in _chebi_dict_of_lists.keys():
                _output_dict_of_lists[_column].append('')

    print_dict_lengths(_uniprot_dict_of_lists, '_uniprot_dict_of_lists')
    print_dict_lengths(_chebi_dict_of_lists, '_chebi_dict_of_lists')
    print_dict_lengths(_output_dict_of_lists, '_output_dict_of_lists')

    # move dictionary of lists into Pandas dataframe at once
    df = pd.DataFrame(_output_dict_of_lists)
    if myoptions.debug: print_df(df)
    write_csv_and_xls(df)
    #if myoptions.outfmt == "csv":
    #    df.to_csv("TPSdownloader.csv", index=False, sep='\t')
    #elif myoptions.outfmt == "xls":
    #    df.to_excel("TPSdownloader.xls", index=False)
    #else:
    #    # maybe more formats would be helpful?
    #    df.to_excel("TPSdownloader.xls", index=False)

if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
