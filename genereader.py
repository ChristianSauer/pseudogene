import io
import re
from bs4 import BeautifulSoup
import pandas as pd
import urllib
import urllib2
from step2 import step2

primers = pd.read_excel("Primer_IGF2.xlsx")
# print primers
primers["Step 1"] = ""
for i, primer in primers.iterrows():
    bad_data = False
    print "Doing {}  of {}".format(i, len(primers))

    if i != 0:
        pass
        # break

    print "Current Primer: {}".format(primer)

    data = {
        "hgsid": "",
        "org": "human",
        "db": "hg19",
        "wp_target": "genome",
        "wp_f": primer["Primer forward"],
        "wp_r": primer["Primer reverse"],
        "Submit": "submit",
        "wp_size": 4000,
        "wp_perfect": 15,
        "wp_good": 15,
        "boolshad.wp_flipReverse": 0
    }

    encoded_data = urllib.urlencode(data)
    content = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/hgPcr", encoded_data)

    soup = BeautifulSoup(content.read())

    try:
        pre = soup.find('pre')
        data = unicode(pre)

        if data.count("<a href") > 1:
            bad_data = True

        data = '\n'.join(data.split('\n')[1:-1])
        # print unicode(data)
    except AttributeError as e:
        print 'No tables found, exiting'

    if bad_data:
        primers["Step 1"][i] = u"Mindestens zwei Matches, Primer unbrauchbar"
        print("Two matches, primer unusable")
        continue

    primers["Step 1"][i] = data

    print "Step 1 done!"

    step2(data, primer["Exon"])

primers.to_excel("result.xlsx")

print "all finished"
