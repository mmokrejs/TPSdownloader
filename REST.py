#! /usr/bin/python3

import io
import os, sys, time
import urllib, requests


def get_get(request):
    _GET = {}
    if request.args != None:
        loc_get = request.args.split("&")
        for stor_key in loc_get:
            tmp = stor_key.split("=")
            if not _GET.has_key(tmp[0]):
                _GET[tmp[0]] = list()
            if len(tmp) > 1:
                _GET[tmp[0]].append(tmp[1])
            else:
                _GET[tmp[0]].append(None)
    return _GET


def get_post(request):
    _dict = request.form
    _POST = {}
    for _key_name in _dict.keys():
        _value = _dict[_key_name]
        parts = _key_name.split('.')
        # initialization 
        _dictionary = _POST
        _key = ""
        _index = 0
        while len(parts) > 0:
            part = parts[0]
            if part[-1:] == "]":
                is_list = 1
                tmp = part.split('[')
                _key = tmp[0]
                _index = int(tmp[1][0:-1])
                if not _dictionary.has_key(_key):
                    _dictionary[_key] = []
                while len(_dictionary[_key]) < _index + 1:
                    _dictionary[_key].append(None)
                if len(parts) > 1:
                    if _dictionary[_key][_index] == None:
                        _dictionary[_key][_index] = {}
                    _dictionary = _dictionary[_key][_index]
            else:
                is_list = 0
                _key = part
                if len(parts) > 1:
                    if not _dictionary.has_key(_key) or _dictionary[_key] == None:
                        _dictionary[_key] = {}
                    _dictionary = _dictionary[_key]
            del(parts[0])
        if is_list == 1:
            _dictionary[_key][_index] = _value
        else:
            _dictionary[_key] = _value
    return _POST


def _q(filename, format='xml', URL='https://www.uniprot.org/uploadlists/'):
    """curl -i -X POST -d "file=@3queries.txt" https://www.uniprot.org/uploadlists/
    curl -i -X POST -H "Content-Type: form-data" -d "file=@3queries.txt" https://www.uniprot.org/uploadlists/
    """

    # In the event you are posting a very large file as a multipart/form-data request, you may want to stream the request. By default, requests does not support this, but there is a separate package which does - requests-toolbelt. You should read the toolbelt’s documentation for more details about how to use it.
    _headers = {'Content_Type': 'form-data', 'User-Agent': 'mmokrejs@gmail.com'}
    if filename and os.path.exists(filename):
        _data = {'files': open(filename, 'rb'), 'format': 'xml', 'from': 'ACC+ID', 'to': 'ACC'}
    else:
        _data = {}

    print("trying URL=%s, data=%s, headers=%s" % (URL, str(_data), str(_headers)))

    # https://docs.python-requests.org/en/master/user/quickstart/#make-a-request
    # https://docs.python-requests.org/en/master/user/advanced/#advanced
    response = requests.post(URL, data=_data, headers=_headers)

    # this does not work with Uniprot
    # POST a Multipart-Encoded File
    # response = requests.post(URL, files={'files': open(filename, 'rb')})

    # <input id="uploadfile" name="file" class="uploadList" type="file"/>
    print("response.headers: %s" % str(response.headers))
    if response.status_code == 200:
        print("response.headers: %s" % str(response.headers))
    elif response.status_code in [ 413, 429, 503 ]:
        # https://developer.mozilla.org/en-US/docs/Web/HTTP/Status#client_error_responses
        time.sleep(int(response.headers.get('retry-after')))
        print("response.headers: %s" % str(response.headers))
        while int(response.headers["Retry-After"]) in [ 413, 429, 503 ]:
            print("response.headers: %s" % str(response.headers))
            time.sleep(int(response.headers.get('retry-after')))
    else:
        sys.stderr.write("Unknown error when writing to %s due to status code %s\n" % (str(URL), response.status_code))

#    except requests.exceptions.Timeout:
#        sys.stderr.write("Timeout when writing to %s due to %s\n" % (str(URL), requests.status_code))
#        return None
#    except:
#        sys.stderr.write("Unknown error when writing to %s due to %s\n" % (str(URL), requests.status_code))
#        return None

    print("response.headers: %s" % str(response.headers))

    # https://urllib3.readthedocs.io/en/latest/reference/urllib3.util.html
    #
    # respect_retry_after_header (bool) – Whether to respect Retry-After header on status codes defined as Retry.RETRY_AFTER_STATUS_CODES or not.
    #
    # get_retry_after(response) - Get the value of Retry-After in seconds.

    print("req: %s" % response)
    print("req.text: %s" % response.text)

    #dir(response)
    #help(response)
    print("Response is:", str(response))

    # id="basket-download-button"

    # from BeautifulSoup import BeautifulSoup
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(response.text, features="lxml")
    links = soup.findAll('a href')
    for link in links:
        print(str(link))

    from pyquery import PyQuery
    pq = PyQuery(response.text)
    _title = pq('title')
    print("Title:", _title)

    # <head><title>yourlist:M2021050472FEB3358BE035486EE75ADE9E917725036DBB9 in UniProtKB</title>
    # https://www.uniprot.org/uniprot/?query=yourlist:M2021050472FEB3358BE035486EE75ADE9E917725036DBB9&format=xml&force=true&sort=yourlist:M2021050472FEB3358BE035486EE75ADE9E917725036DBB9&compress=yes

    # if response.ok:
    #     _res = response.read().decode("utf8") # AttributeError: 'Response' object has no attribute 'read'
    #     dir(_res)
    #     print(str(_res))
    # else:
    #     _res = None

    return response.text


def main():
    _h = _q('3queries.txt')
    print(str(_h))


if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
