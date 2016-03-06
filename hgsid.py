from bs4 import BeautifulSoup
import requests

__author__ = 'Christian'


def get_hgsid():
    with open("hgsid.txt", "r") as hgsid_source:
        return hgsid_source.read()


# def get_hgsid(url):
#    hgsid = hgsid_from_response(requests.get(url))
#    return hgsid


def hgsid_from_response(r):
    soup = BeautifulSoup(r.text)
    # print r.text
    hgsid = set([i.get('value') for i in soup('input', dict(name='hgsid'))])
    assert len(hgsid) == 1
    hgsid = list(hgsid)[0]
    return hgsid