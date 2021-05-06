#! /usr/bin/python3

import io
import os, sys
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
    _headers = {'Content_Type': 'form-data', 'User-Agent': 'mmokrejs@gmail.com'}
    if filename and os.path.exists(filename):
        _data = {"file": open(filename, 'rb')}
    else:
        _data = {}

    print("trying URL=%s, data=%s, headers=%s" % (URL, str(_data), str(_headers)))
    try:
        _req = requests.post(URL, data=_data, headers=_headers, timeout=60)
    except requests.exceptions.Timeout:
        sys.stderr.write("Timeout when writing to %s\n" % str(URL))
        return None
    except:
        sys.stderr.write("Unknown error when writing to %s\n" % str(URL))
        return None

    print("req: %s" % _req)

    #dir(_req)
    #help(_req)
    print("Response is:", str(_req))

    if _req.ok:
        _res = _req.read().decode("utf8")
        dir(_res)
        print(str(_res))
    else:
        _res = None

    return _res


def main():
    _h = _q('3queries.txt')
    print(str(_h))


if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
