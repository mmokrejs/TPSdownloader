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


Mass-download from Uniprot

Currently, I manually downloaded the entries from UniProt as `xml.gz` files.
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
