# TPSdownloader
Fetch TPS cyclases information from various sources

This utility can parse into pythonic structures data from XLS table,
see https://docs.google.com/spreadsheets/d/1ymHBHuy5Dwzi_RvK233tHtfXWKnNTWTZtwTSV1ALVHE/edit?skip_itp2_check=true#gid=0

Then the utility parses pre-downloaded XML data from a local cache. The files
can be typically single-entry uncompressed XML files or compressed multi-entry
XML files from UniProt.

Then, this utility fetches XML entries for UniProt accessions specified on
a command-line or in a single-column TXT file. If the entries point to ChEBI,
it fetches them too. It does not download some entries already placed on a blacklist,
like water, some ions, etc. If they get fetched from ChEBI, the get cached in a local
file but then, nevertheless, are not used further.

Uniprot has obsoleted lots of duplicated protein in 2015 originating
from closely related proteomes. Therefore, some proteins has been
withdrawn from Uniprot. Such Uniprot IDs not fetched&processed are
listed at the end of the run as obsoleted protein IDs. For more info
please refer to https://www.uniprot.org/help/proteome_redundancy


*Requirements*

pandas
natsort
xlsxwriter
conda
rdkit


*Usage*

The program can be run in several ways. It always outputs several files, e.g.:
```
Info: Wrote TPSdownloader_20210717172721.xlsx file.
Info: Wrote TPSdownloader_20210717172721.csv file.
Info: Wrote TPSdownloader_product_ChEBI_IDs_20210717172721.csv file.
Info: Wrote TPSdownloader_EC_numbers_20210717172721.csv file.
Info: Wrote TPSdownloader_Rhea_IDs_20210717172721.csv file.
```

It also outputs a list of ChEBI IDs (substrates, products, cofactors, junk), a list product ChEBI IDs, a list of associated EC and Rhea numbers.

1.

It can fetch IDs from network, parse them.

`python TPSdownloader.py --uniprot-id-file=Uniprot-Uniparc_domain-hits-to-terpene_synth_domains.txt`


2.

It can also mix the fetched data with XLSX data already annotated
(while this does NOT preserve all the manual data).

`python TPSdownloader.py --uniprot-id-file=Uniprot-Uniparc_domain-hits-to-terpene_synth_domains.txt --xls-storage=TPS-database_20210717.xlsx`

The above approach will detect only 72 new proteins in the parsed data which we still lack in a manually curated table, of those 7 are annotated in Uniprot with ChEBI/Rhea IDs or EC numbers or description of the reaction at least. It outputs 25161 rows in the main table, multiple lines per single Uniprot ID, while preeserving such multiplicated input rows (if substrates and products were different).


3.

It can also parse already fetched XML files, parse it, process a list
of Uniprot IDs to be transferred from the parsed data into resulting files.
Optionally, it can take as an input a list of Uniprot IDs we have already
curated and although the proteins will appear in the output table, their
proteins will NOT appear in the list of distinct proteins available in parsed
data (see XLSX sheet `New proteins to be curated`). The new two command are
practically doing the same. The first-one grabs Uniprot IDs from such column
in the XLSX file, the second using a simple TXT file with the Uniprot IDs.

`python TPSdownloader.py --uniprot-ids-from-file=TPS-database_20210717.xlsx --already-curated-id-from-file=Uniprot_IDs_manually_curated.txt`
`python TPSdownloader.py --uniprot-id-file=Uniprot-Uniparc_domain-hits-to-terpene_synth_domains.txt --already-curated-id-from-file=Uniprot_IDs_manually_curated.txt`


*Bugs*

For an unknown reason the command

`python TPSdownloader.py --uniprot-ids-from-file=TPSdownloader_20210717135125.xlsx --already-curated-id-from-file=Uniprot_IDs_manually_curated.txt`

results in 160 new proteins in the parsed data which we still lack in a manually curated table, of those 12 are annotated in Uniprot with ChEBI/Rhea IDs or EC numbers or description of the reaction at least. It outputs 24799 rows in the main table, only.

For example, the entry `H8ZM70` appears on a single row although it should be split across 4 rows as it catalyzes 4 reactions:
```
[['(2E,6E,10E)-geranylgeranyl diphosphate = (+)-copalyl diphosphate',
  '(+)-copalyl diphosphate = abieta-7,13-diene + diphosphate',
  '(+)-copalyl diphosphate = diphosphate + neoabietadiene',
  '(+)-copalyl diphosphate = abieta-8(14),12-diene + diphosphate']]
```

In the manually curated table we have all these 4 reactions recorded.

Improving parsing the `--xls-storage` contents would need a lot of changes.



*Mass-download from Uniprot*

Currently, I manually downloaded the entries from UniProt as `xml.gz` files,
the parser was able to use *.xml.gz files but after switching from
`elementtree.parse()` to to a iterable parser `elementtree.iterparse()` we
lost the ability to act of decompressed gzip stream as the decompressed file
object lacks some properties. So the files must be uncompressed as of now.

`split --lines 4000 --numeric-suffixes=0 Uniprot_ID_od_Terezy.tsv queries.batch.`


Answer from UniProt

```
Dear Martin,

The correct curl syntax is

curl -X POST
      -L
      "https://www.uniprot.org/uploadlists/"
      -F "file=@3queries.txt"                   #note the file name
      -F "format=xml"
      -F "from=ACC+ID"
      -F "to=ACC"

Different examples for ID mapping are available at.

https://www.uniprot.org/help/api%5Fidmapping

Of which batch retrieval is a sub usecase.

To note this API sends a redirect (curl -L) python likes
to take exception to this so be prepared to catch the Exception and
follow the redirect.

Regards,
Jerven


On 07/05/2021 09:18, Elisabeth Gasteiger via RT wrote:
>
> Fri May 07 09:18:29 2021: Request 209412 was acted upon.
> Transaction: Given to jbolleman (Jerven Bolleman) by gasteig
>         Queue: UniProt
>       Subject: [uuw] Is web REST api really working?
>         Owner: jbolleman
>    Requestors: mmokrejs@fold.natur.cuni.cz
>        Status: new
>   Ticket <URL: https://rt.isb-sib.ch/Ticket/Display.html?id=209412 >
>
>
> This transaction appears to have no content
>
>
> Initial request :
> ----------------------------------------------------------------
>
> Hi,
>    I would like to fecth about 27k entries from Uniprot using remote API in XML format. Seems https://www.uniprot.org/help/api_batch_retrieval and https://www.uniprot.org/help/api_batch_retrieval is supposed to work, although you are very sparse in explaining what are the key bits to be eventually re-implememted by me in say python. After many hours I am giving, because not even these work using curl standalone utility:
>
> curl -X POST -d "name=@3queries.txt" https://www.uniprot.org/uploadlists/
> curl -X POST -H "Content-Type: form-data" -d "@3queries.txt" https://www.uniprot.org/uploadlists/
> curl -i -X POST -H "Content-Type: form-data" -d "@3queries.txt" https://www.uniprot.org/uploadlists/
> curl -X POST -F "file=@3queries.txt" https://www.uniprot.org/uploadlists/
>
>
> Here is a hacking code in python, which I conclude is not at fault at the moment
>
> def _q(filename, format='xml', URL='https://www.uniprot.org/uploadlists/'):
>      """curl -i -X POST -d "name=@3queries.txt" https://www.uniprot.org/uploadlists/
>      curl -i-X POST -H "Content-Type: form-data" -d "@3queries.txt" https://www.uniprot.org/uploadlists/
>      """
>
>      _headers = {'Content-Type': 'form-data', 'User-Agent': 'mmokrejs@gmail.com'}
>      if filename and os.path.exists(filename):
>          _data = {"file": open(filename, 'rb'), 'format': 'xml', 'from': 'ACC+ID', 'to': 'ACC'}
>      else:
>          _data = {}
>
>      print("trying URL=%s, data=%s, headers=%s" % (URL, str(_data), str(_headers)))
>      try:
>          _req = requests.post(URL, data=_data, headers=_headers, timeout=60)
>      except requests.exceptions.Timeout:
>          sys.stderr.write("Timeout when writing to %s due to %s\n" % (str(URL), requests.status_code))
>          return None
>      except:
>          sys.stderr.write("Unknown error when writing to %s due to %s\n" % (str(URL), requests.status_code))
>          return None
>
>      print("_req.headers: %s" % str(_req.headers))
>
>      while 'Retry-After' in str(_req.headers):
>          print("sleeping for 2 secs")
>          sleep(2)
>      print("req: %s" % _req)
>      print("req.text: %s" % _req.text)
>      print("Response is:", str(_req))
>      return _req.text
>
>
> def main():
>      _h = _q('3queries.txt') # file with 3 Uniprot ID's, line by line
>      print(str(_h))
>
>
> if __name__ == "__main__":
>      main()
>
>
> Thank you for your help,
> Martin
>
> Name: Martin Mokrejs
> Referred from: https://www.uniprot.org/help/about
> Browser: Mozilla/5.0 (X11; Linux x86_64; rv:87.0) Gecko/20100101 Firefox/87.0
> IP address: 89.24.41.226
>
>
>
>
> Complete Ticket History :
> ----------------------------------------------------------------
>

-- 

	*Jerven Tjalling Bolleman*
Principal Software Developer
*SIB | Swiss Institute of Bioinformatics*
1, rue Michel Servet - CH 1211 Geneva 4 - Switzerland
t +41 22 379 58 85
Jerven.Bolleman@sib.swiss - www.sib.swiss
```
