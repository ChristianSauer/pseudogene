import os
import urllib
import urllib2
from bs4 import BeautifulSoup
import re
import io
import pandas as pd
__author__ = 'Christian'


def step2(gene_sequence, exon, steps_to_excel=True, local_path=""):
    print("Starting Step 2")
    # Step 2
    data = {
        "hgsid": "",
        "changeInfo": "",
        "org": "Human",
        "db": "hg19",
        "type": "BLAT's guess",
        "sort": "query,score",
        "output": "hyperlink",
        "userSeq": gene_sequence,
        "Submit": "submit",
    }

    encoded_data = urllib.urlencode(data)
    content = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/hgBlat", encoded_data)

    soup = BeautifulSoup(content.read())
    # and again in a pre tag
    try:
        pre = soup.find('pre')
        data = unicode(pre)
        data = data.replace("<pre>", "").replace("</pre>", "")
        data = data.replace("<a href", "<a_href")
        data = "Browser " + data
        data = re.sub("--+", "", data)
        data_sane = pd.read_csv(io.StringIO(unicode(data)), sep=r"\s+")
        data_sane.reset_index()

        data_sane.rename(columns={'ACTIONS': 'Details'}, inplace=True)

        for col_name in ["Browser", "Details"]:
            def to_usable_url(url_str):
                assert isinstance(url_str, str)
                url_str = url_str.replace('<a_href="..', "http://genome.ucsc.edu/")
                right_end = url_str.find('">')
                url_str = url_str[:right_end]

                if len(url_str) > 255:
                    # stupid excel limitation of 255 chars in single formula
                    url_str = "SPLITME " + url_str

                return url_str

            data_sane[col_name] = data_sane[col_name].apply(to_usable_url)
            data_sane["Exon"] = exon

        filename = "{}.xlsx".format(exon)
        local_path = os.path.join(local_path, filename)
        print("Step 2 done for: " + exon)

        if steps_to_excel:
            data_sane.to_excel(local_path)
        return data_sane
    except AttributeError:
        print 'No tables found, exiting'

