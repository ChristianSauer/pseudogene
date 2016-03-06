# coding=utf-8
from StringIO import StringIO
import json
import os
import requests
from step2 import step2
from bs4 import BeautifulSoup
import pandas as pd
import multiprocessing

__author__ = 'Christian Sauer'

pd.options.mode.chained_assignment = None


def prepare_input():
    lines = open(source_file).readlines()
    lines = lines[config["skip_first_n_lines"]:]

    pos = -1
    for line in lines:
        pos += 1
        if line.startswith(config["skip_lines_after"]):
            lines = lines[pos + 1:]
            break
    else:
        print ("no duplicate tracks found")

    column_names = [x.strip() for x in config["column_names"].split(",")]

    io_lines = StringIO(u"\n".join(lines))
    all_primers = pd.read_table(io_lines, names=column_names, header=None)

    print "Read %s primers" % len(all_primers)

    return all_primers


def find_genes():
    found_genes_df = primers[primers['names'].str.contains(search_gene)]
    if len(found_genes_df) == 0:
        raise IOError("No genes selected")

    print("Genes found: \n" + "\n".join(found_genes_df['names'].unique()))
    print "found %s genes" % len(found_genes_df)

    return found_genes_df


def collect_results():
    # assumption: first gene is always the main - match.
    # We are trying to find duplicates (pseudogenes) in the rows after that
    # a pseudogene has a high identity and similar span

    filtered_dfs = []
    for df in results:
        print "Checking Exon:{}".format(df["Exon"][0])
        reference_row = df.iloc[0]
        df["span_similarity"] = df.apply(lambda x: (100 * x["SPAN"]) / reference_row["SPAN"], axis=1)
        df["identity_number"] = df.apply(lambda x: float(x["IDENTITY"].strip().replace("%", "")), axis=1)
        identity_mask = df["identity_number"] >= config["min_identity"]
        min_similarity_mask = df["span_similarity"] >= config["min_span_similarity"]
        max_similarity_mask = df["span_similarity"] < config["max_span_similarity"]

        combined_mask = identity_mask & min_similarity_mask & max_similarity_mask

        # always ignore first row
        combined_mask[0] = False

        print "{} potential pseudogenes found.".format(sum(combined_mask))
        matching_df = df[combined_mask]
        matching_df.drop('identity_number', axis=1, inplace=True)
        filtered_dfs.append(matching_df)

    master_df = filtered_dfs.pop(0)
    for df in filtered_dfs:
        master_df = master_df.append(df)

    master_df["config_min_identity"] = config["min_identity"]
    master_df["config_min_span_similarity"] = config["min_span_similarity"]
    master_df["config_max_span_similarity"] = config["max_span_similarity"]
    xlsx_name = "Pseudognes for {}.xlsx".format(search_gene)
    local_path = os.path.join(folder_name, xlsx_name)
    print "Writing pseudogenes into: " + local_path
    master_df.to_excel(local_path)


def prepare_reduced_df():
    def create_pos19_str(x):
        return x["chromosome"] + ":" + str(x["start [hg19]"]) + "-" + str(x["end [hg19]"])

    reduced_df["position [hg19]"] = reduced_df.apply(create_pos19_str, axis=1)
    reduced_df["gene name"] = search_gene
    reduced_df["sequence length"] = reduced_df["end [hg19]"] - reduced_df["start [hg19]"]
    reduced_df["genomic sequence [hg19]"] = ""


def get_data(args):
    current_row, row, hg_sid, cfg, folder = args

    dna_pos = row["position [hg19]"]

    try:
        data = {"hgsid": hg_sid,
                "table": "",
                "o": "16370229",
                "getDnaPos": dna_pos,
                "hgSeq.cdsExon": "1",
                "hgSeq.padding5": "0",
                "hgSeq.padding3": "0",
                "hgSeq.casing": "upper",
                "boolshad.hgSeq.maskRepeats": "0",
                "hgSeq.repMasking": "lower",
                "boolshad.hgSeq.revComp": "0",
                "submit": "get DNA",
                }

        r = requests.get("http://genome.ucsc.edu/cgi-bin/hgc?g=htcGetDna2", params=data, allow_redirects=True)

        soup = BeautifulSoup(r.content)

        try:
            pre = soup.find(u'pre')
            data = unicode(pre.getText())
            data = data.splitlines()
            data = data[2:]
            data = u"".join(data)
            print(data)
        except AttributeError:
            print 'No tables found, exiting'

        name = "{} {} {}-{}".format(row["gene name"], current_row, row["start [hg19]"], row["end [hg19]"])

        results_df = step2(data, name, cfg["detailed_results_to_excel"], folder)

        return results_df, data, current_row
    except Exception as ex:
        print(ex)
        print(
            "Could not execute getDNA for {}, HGSID stale? Get a new one and paste it into config.json".format(
                dna_pos))


def get_and_append_dna():
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() * 2)

    results_for_pseudogenes = []

    def get_result(args):
        result, data, current_row = args
        # This is called whenever foo_pool(i) returns a result.
        # result_list is modified only by the main process, not the pool workers.
        reduced_df.set_value(current_row, "genomic sequence [hg19]", data)
        results_for_pseudogenes.append(result)

    # result = [get_data((i, row)) for i, row in reduced_df.iterrows()]
    result = pool.map(get_data, ((i, row, hgsid, config, folder_name) for i, row in reduced_df.iterrows()))
    _ = [get_result(x) for x in result if x is not None]
    pool.close()
    pool.join()

    return results_for_pseudogenes


def get_folder_name(gene_name):
    import shutil
    local_folder = os.path.join(os.getcwd(), gene_name)
    if os.path.isdir(local_folder):
        msg = "Warning Folder {} already exists. Deleting all content and continue? y or n".format(local_folder)
        if raw_input(msg).lower().strip() == "n":
            exit()
        else:

            shutil.rmtree(local_folder)

    os.makedirs(local_folder)
    return local_folder


if __name__ == "__main__":
    with open('config.json') as data_file:
        config = json.load(data_file)

    hgsid = config["hgsid"]
    print ("Using hgsid: " + hgsid)

    source_file = config["source"]
    print "Source File: " + source_file

    primers = prepare_input()

    print "Please enter a gene name"

    search_gene = raw_input('--> ')
    search_gene = search_gene.strip().upper()

    print "searching for: %s" % search_gene

    reduced_df = find_genes()

    if raw_input("Continue? y or n").lower().strip() == "n":
        exit()

    folder_name = get_folder_name(search_gene)

    prepare_reduced_df()
    results = get_and_append_dna()

    print "Step 3: collecting results into output table"
    collect_results()

    file_name = "DNA FOR {}.xlsx".format(search_gene)
    local_result_path = os.path.join(folder_name, file_name)
    reduced_df.to_excel(local_result_path)
    print("Written to: {}".format(local_result_path))
    print("Done")
