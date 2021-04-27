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

version = "20210422"
myparser = OptionParser(version="%s version %s" % ('%prog', version))
myparser.add_option("--uniprot-id-file", action="store", type="string", dest="uniprot_idfile", default='-',
    help="List of Uniprot IDs in a file or STDIN, is mutually exclusive with --chebi-id-file.")
myparser.add_option("--chebi-id-file", action="store", type="string", dest="chebi_idfile", default='-',
    help="List of ChEBI IDs in a file or STDIN, is mutually exclusive with --uniprot-id-file.")
myparser.add_option("--uniprot-id", action="store", type="string", dest="uniprot_id", default=None,
    help="List of Uniprot IDs, is mutually exclusive with --chebi-id.")
myparser.add_option("--chebi-id", action="store", type="string", dest="chebi_id", default=None,
    help="List of ChEBI IDs, is mutually exclusive with --uniprot-id.")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug to some value")
(myoptions, myargs) = myparser.parse_args()


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


def parse_chebi_xml(filename):
    """Keeps data in lists so we can capture multiple values eventually, do not know
    what to expect at the moment."""

    etree=ET.parse(filename)
    root=etree.getroot()

    _chebi_id = []
    _definition = []
    _names = []
    _formula = []
    _smiles = []
    _inchi = []

    for elem in root:
        if myoptions.debug: print("Level 0: ", elem.tag, ' ', elem.attrib, ' ', elem.text)
        for child in elem:
            if myoptions.debug > 1: print("Items: ", str(child.items()))

            if myoptions.debug: print("Level 1: tag:", child.tag, 'attrib:', child.attrib, 'text:', child.text)
            if child.tag == 'ID':
                if myoptions.debug: print("L1: child.attrib: ", child.attrib, "child.tag: ", child.tag)
                if not _chebi_id:
                    _chebi_id = [child.text]
                else:
                    raise(ValueError)
            if child.tag == 'NAME':
                if not _names:
                    _names = [child.text]
                else:
                    _names += [child.text]
            if child.tag == 'DEFINITION':
                if not _definition:
                    _definition = [child.text]
                else:
                    _definition += [child.text]
            if child.tag == 'FORMULA':
                if not _formula:
                    _formula = [child.text]
            if child.tag == 'SYNONYM':
                if not _names:
                    _names = [child.text]
                else:
                    _names += [child.text]
            if child.tag == 'SMILES':
                if not _smiles:
                    _smiles = [child.text]
            if child.tag == 'INCHI':
                if not _inchi:
                    _inchi = [child.text]

    print("Info: IDs: %s, names: %s, definition: %s, formula: %s, smiles: %s, inchi: %s" % (str(_chebi_id), str(_names), str(_definition), str(_formula), str(_smiles), str(_inchi)))
    return(_chebi_id, _names, _definition, _formula, _smiles, _inchi)


def parse_uniprot_xml(filename, terpenes):
    etree=ET.parse(filename)
    root=etree.getroot()

    _accessions = []
    _uniprot_name = None
    _protein_names = [] # sometimes UniProt has unset protein_names and also feature_description
    _feature_descriptions = []
    _chebi_ids = []
    _sequence = None
    _organism = None
    _lineage = None

    # blacklist of IDs of reaction substrates or non-cyclic terpenes, etc., but including β-farnesene CHEBI:10418
    non_terpene_chebi_ids = ['CHEBI:15377', 'CHEBI:33019', 'CHEBI:58057', 'CHEBI:128769', 'CHEBI:175763', 'CHEBI:57533', 'CHEBI:10418', 'CHEBI:10280', 'CHEBI:175763', 'CHEBI:64283', 'CHEBI:58756', 'CHEBI:15441', 'CHEBI:58622', 'CHEBI:58553', 'CHEBI:57665', 'CHEBI:58635', 'CHEBI:15440', 'CHEBI:138223', 'CHEBI:128769', 'CHEBI:57907', 'CHEBI:64801', 'CHEBI:61984', 'CHEBI:58206', 'CHEBI:138167', 'CHEBI:15347', 'CHEBI:162247', 'CHEBI:17221', 'CHEBI:60374', 'CHEBI:61746', 'CHEBI:98']
    # CHEBI:33019 - diphosphate(3−)
    # CHEBI:58057 - geranyl diphosphate(3−)
    # CHEBI:128769 - isopentenyl diphosphate(3−)
    # CHEBI:175763 - 2-trans,6-trans-farnesyl diphosphate(3−)
    # CHEBI:57533 - geranylgeranyl diphosphate(3−)
    # CHEBI:15377 - water
    
    for elem in root:
        if myoptions.debug: print("Level 0: ", elem.tag, ' ', elem.attrib, ' ', elem.text)
        for child in elem:
            if myoptions.debug > 1: print("Items: ", str(child.items()))
    
            if myoptions.debug: print("Level 1: tag:", child.tag, 'attrib:', child.attrib, 'text:', child.text)
            if child.tag == '{http://uniprot.org/uniprot}accession':
                if myoptions.debug: print("L1: child.attrib: ", child.attrib, "child.tag: ", child.tag)
                if not _accessions:
                    _accessions = [child.text] # "Q6XDB5"
                else:
                    _accessions += [child.text] # "C0PT91"
            elif child.tag == '{http://uniprot.org/uniprot}name':
                _uniprot_name = child.text # "TPSD2_PICSI"
            elif child.tag == '{http://uniprot.org/uniprot}sequence':
                _sequence = child.text
            if child.tag == '{http://uniprot.org/uniprot}feature':
                if 'type' in child.attrib.keys() and 'description' in child.attrib.keys() and child.attrib['type'] == 'chain':
                    _feature_descriptions += [child.attrib['description']]
                    
            for subchild in child:
                if myoptions.debug: print("Level 2: ", subchild.tag, ' ', subchild.attrib, ' ', subchild.text)
                _chebi_ids_local = []
                if child.tag == '{http://uniprot.org/uniprot}organism':
                    if subchild.tag == '{http://uniprot.org/uniprot}name' and subchild.attrib['type'] == 'scientific':
                        _organism = subchild.text
                for sschild in subchild:
                    tag = {}
                    if myoptions.debug: print("Level 3: ", sschild.tag, ' ', sschild.attrib, ' ', sschild.text)
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
                            if not _chebi_ids_local:
                                _chebi_ids_local = [sschild.attrib['id']]
                            else:
                                _chebi_ids_local += [sschild.attrib['id']]

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
                    _chebi_ids += [_chebi_ids_local]
                    _chebi_ids_local = []

    if myoptions.debug:
        for _i in range(0,len(_chebi_ids)-1):
            if len(_chebi_ids) != len(_protein_names):
                print("Warning: Number of ChEBI entries does not match number of alternative protein names : _chebi_ids=%s _protein_names=%s" % (str(_chebi_ids), str(_protein_names)))
                print("Info: accessions: %s, chebi_ids: %s, _protein_names: %s, description: %s, organism: %s, lineage: %s, sequence: %s" % (str(_accessions), str(_chebi_ids[_i]), str(_protein_names), str(_feature_descriptions), str(_organism), str(_lineage), str(_sequence)))
            else:
                print("Info: accessions: %s, chebi_ids: %s, _protein_names: %s, description: %s, organism: %s, lineage: %s, sequence: %s" % (str(_accessions), str(_chebi_ids[_i]), str(_protein_names[_i]), str(_feature_descriptions), str(_organism), str(_lineage), str(_sequence)))

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
    return(_accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence)


def create_cache(cachedir=".TPSdownloader_cache"):
    "make local caching directory for input XML files"

    if not os.path.exists(cachedir):
        os.mkdir(cachedir)
    for subdir in ('uniprot', 'chebi'):
        if not os.path.exists(cachedir + os.path.sep + subdir):
            os.mkdir(cachedir + os.path.sep + subdir)


def downloader(cmdline):
    print(cmdline)
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
            _cmdline = "curl --no-progress-meter -o " + _filename + " " + url + myid
            if dbname == 'uniprot':
                _cmdline += ".xml"
            downloader(_cmdline)
            # could also prettyprint the XML files using
            # xml_pp -i.bak _filename
        else:
            sys.stderr.write("Info: File %s already exists\n" % _filename)
    else:
        sys.stderr.write("Error: Directory %s does not exist\n" % _cachedir)


def download_uniprot(myid, path=".TPSdownloader_cache" + os.path.sep + 'uniprot' + os.path.sep):
    """Download a page like https://www.uniprot.org/uniprot/A0A2K9RFZ2.xml

    Some entries have multiple Accessions, like

    <accession>Q6XDB5</accession>
    <accession>C0PT91</accession>
    <name>TPSD2_PICSI</name>

    """

    downloader_wrapper(myid, 'uniprot', ".TPSdownloader_cache" + os.path.sep, "https://www.uniprot.org/uniprot/")


def download_chebi(myid, path=".TPSdownloader_cache" + os.path.sep + 'chebi' + os.path.sep):
    """Download a page like https://www.ebi.ac.uk/chebi/saveStructure.do?xml=true&chebiId=58622

    ChEBI also provides SQL dumps, probably an overkill for our
    purpose.
    """

    downloader_wrapper(myid, 'chebi', ".TPSdownloader_cache" + os.path.sep, "https://www.ebi.ac.uk/chebi/saveStructure.do?xml=true&chebiId=")


def import_existing_csv(filename="TPS-database_20210420.xlsx",sheetname="Sheet1"):
    """Import existing table and add CHEBI Id column, InChI and InChIKey
    if missing in our table and eventually already recorded in Uniprot.
    """
    
    df = pd.read_excel(filename, sheetname, index_col=None, na_values=["NA"])
    return df


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


def write_csv(df):
    """Write CSV output, per one Uniprot ID A0A2K9RFZ2 output even
    multiple lines if there are multiple reactions catalyzed

    https://www.uniprot.org/uniprot/A0A2K9RFZ2
    https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:58622
    """
    
    _datetime = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    df.to_excel("TPSdownloader_" + _datetime + ".xlsx", sheet_name="Sheet1")


def fetch_ids_from_xlsx(terpenes):
    _all_uniprot_ids = set()
    _all_chebi_ids = set()
    df = import_existing_csv() # parse UniProt IDs from filename="TPS-database_20210420.xlsx",sheetname="Sheet1"
    for i in df['Uniprot ID']:
        if i is np.NaN:
            _i = None
        elif i is None:
            _i = None
        else:
            _i = i

        if i and isinstance(i, str):
            _i = i.strip()

        if _i is np.NaN:
            pass
        elif _i is None:
            pass
        elif _i == 'nan':
            pass
        elif str(_i) == 'nan':
            pass
        if _i:
            download_uniprot(_i)
            _filename = ".TPSdownloader_cache" + os.path.sep + 'uniprot' + os.path.sep + _i + '.xml'
            if os.path.exists(_filename) and os.path.getsize(_filename):
                (_accessions, _chebi_ids, _protein_names, _feature_descriptions, _organism, _lineage, _sequence) = parse_uniprot_xml(_filename, terpenes)
                for _uniprot_id in _accessions:
                    _all_uniprot_ids.update([_uniprot_id])
                for _chebi_id in _chebi_ids:
                    _all_chebi_ids.update(set(_chebi_id))
        if myoptions.debug or True: print("Debug: _all_uniprot_ids=%s" % str(_all_uniprot_ids))
        if myoptions.debug or True: print("Debug: _all_chebi_ids=%s" % str(_all_chebi_ids))
    return (_all_uniprot_ids, _all_chebi_ids)


_r1 = re.compile(r'C[0-9]+')

def classify_terpene(formula):
    _match = _r1.search(formula)
    if _match:
        _carbon_count = int(formula[_match.span()[0]+1:_match.span()[1]])
    else:
        _carbon_count = None
    if _carbon_count > 12 and _carbon_count < 17:
        _terpene_type = 'sesq'
    elif _carbon_count > 8 and _carbon_count < 12:
        _terpene_type = 'mono'
    elif _carbon_count > 18 and _carbon_count < 22:
        _terpene_type = 'di'
    elif _carbon_count > 28 and _carbon_count < 32:
        _terpene_type = 'tri'
    elif _carbon_count > 38 and _carbon_count < 42:
        _terpene_type = 'tetra'
    elif _carbon_count > 23 and _carbon_count < 27:
        _terpene_type = 'sest'
    else:
        _terpene_type = 'unexpected'
    return _terpene_type


def main():
    create_cache()
    _terpenes = parse_known_terpenes()
    if myoptions.uniprot_id:
        download_uniprot(myoptions.uniprot_id)
        _all_uniprot_ids = set()
        _all_chebi_ids = set()
    else:
        (_all_uniprot_ids, _all_chebi_ids) = fetch_ids_from_xlsx(_terpenes)

    for _chebi_id in _all_chebi_ids:
        download_chebi(_chebi_id)
        _filename = ".TPSdownloader_cache" + os.path.sep + 'chebi' + os.path.sep + _chebi_id + '.xml'
        if os.path.exists(_filename) and os.path.getsize(_filename):
            _chebi_id2, _names, _definition, _formula, _smiles, _inchi = parse_chebi_xml(_filename)
            _terpene_type = classify_terpene(_formula[0])

        #for _i in range(0,len(_chebi_ids)-1):
        #    if len(_chebi_ids) != len(_protein_names):
        #        print("Warning: Number of ChEBI entries does not match number of alternative protein names : _chebi_ids=%s _protein_names=%s" % (str(_chebi_ids), str(_protein_names)))
        #        print("Info: accessions: %s, chebi_ids: %s, _protein_names: %s, description: %s, organism: %s, lineage: %s, sequence: %s" % (str(_accessions), str(_chebi_ids[_i]), str(_protein_names), str(_feature_descriptions), str(_organism), str(_lineage), str(_sequence)))
        #    else:
        #        print("Info: accessions: %s, chebi_ids: %s, _protein_names: %s, description: %s, organism: %s, lineage: %s, sequence: %s" % (str(_accessions), str(_chebi_ids[_i]), str(_protein_names[_i]), str(_feature_descriptions), str(_organism), str(_lineage), str(_sequence)))

    # move dictionary of lists into Pandas dataframe at once
    # dfObj = pd.DataFrame({'Uniprot ID': _accessions, 'Name': _protein_names, 'Amino acid sequence': _sequence, 'Species': _organism, 'Kingdom': _lineage, 'Terpene type': _terpene_type, 'Substrate (including stereochemistry)': None, 'Cofactors': None, 'Name of intermediate': None, 'SMILES of intermediate': None, 'Product': None, 'SMILES of product (including stereochemistry)': None, 'CheEBI ID': _chebi_ids, 'Notes': None, 'Publication (URL)': None})

if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
