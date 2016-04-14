# coding=utf-8
from StringIO import StringIO
import json
import logging
import os
import datetime
import shutil
import requests
from step2 import step2
from bs4 import BeautifulSoup
import pandas as pd
import multiprocessing

__author__ = 'Christian Sauer'

pd.options.mode.chained_assignment = None
logging.basicConfig(level=logging.DEBUG, filename="log.txt", filemode="a+",
                    format="%(asctime)-15s %(levelname)-8s %(message)s")

logger = logging.getLogger(__name__)


class PseudoGeneFinder(object):
    def __init__(self, config=None, genes=None):

        if config is None:
            config_path = 'config.json'
        else:
            config_path = config

        with open(config_path) as data_file:
            self.config = json.load(data_file)

        self.hgsid = self.config["hgsid"]
        self.source_file = self.config["source"]
        self.primers = None
        self.reduced_df = None
        self.search_gene = None
        self.folder_name = None
        self.results = None
        self.genes = genes
        self._unattended = genes is not None

    def _prepare_input(self):
        lines = open(self.source_file).readlines()
        lines = lines[self.config["skip_first_n_lines"]:]

        pos = -1
        for line in lines:
            pos += 1
            if line.startswith(self.config["skip_lines_after"]):
                lines = lines[pos + 1:]
                break
        else:
            logger.info("no duplicate tracks found")

        column_names = [x.strip() for x in self.config["column_names"].split(",")]

        io_lines = StringIO(u"\n".join(lines))
        all_primers = pd.read_table(io_lines, names=column_names, header=None)

        logger.info("Read %s primers" % len(all_primers))

        return all_primers

    def _find_genes(self):
        found_genes_df = self.primers[self.primers['names'].str.contains(self.search_gene)]
        if len(found_genes_df) == 0:
            raise IOError("No genes selected")

        logger.info("Genes found: \n" + "\n".join(found_genes_df['names'].unique()))
        logger.info("found %s genes" % len(found_genes_df))

        return found_genes_df

    def _collect_results(self):
        # assumption: first gene is always the main - match.
        # We are trying to find duplicates (pseudogenes) in the rows after that
        # a pseudogene has a high identity and similar span

        filtered_dfs = []
        for df in self.results:
            logger.info("Checking Exon:{}".format(df["Exon"][0]))
            reference_row = df.iloc[0]
            df["span_similarity"] = df.apply(lambda x: (100 * x["SPAN"]) / reference_row["SPAN"], axis=1)
            df["identity_number"] = df.apply(lambda x: float(x["IDENTITY"].strip().replace("%", "")), axis=1)
            identity_mask = df["identity_number"] >= self.config["min_identity"]
            min_similarity_mask = df["span_similarity"] >= self.config["min_span_similarity"]
            max_similarity_mask = df["span_similarity"] < self.config["max_span_similarity"]

            combined_mask = identity_mask & min_similarity_mask & max_similarity_mask

            # always ignore first row
            combined_mask[0] = False

            logger.info("{} potential pseudogenes found.".format(sum(combined_mask)))
            matching_df = df[combined_mask]
            matching_df.drop('identity_number', axis=1, inplace=True)
            filtered_dfs.append(matching_df)

        master_df = filtered_dfs.pop(0)
        for df in filtered_dfs:
            master_df = master_df.append(df)

        # always show at least one line with the collecting information
        if len(master_df) == 0:
            master_df = self._append_no_pseudogenes_info(master_df)

        master_df["config_min_identity"] = self.config["min_identity"]
        master_df["config_min_span_similarity"] = self.config["min_span_similarity"]
        master_df["config_max_span_similarity"] = self.config["max_span_similarity"]

        xlsx_name = "Pseudognes for {}.xlsx".format(self.search_gene)
        local_path = os.path.join(self.folder_name, xlsx_name)
        logger.info("Writing pseudogenes into: " + local_path)
        master_df.to_excel(local_path)

    def _append_no_pseudogenes_info(self, master_df):
        dummy_data = {"Browser": "", "Details": "No Pseudogenes found", "QUERY": "", "SCORE": "", "START": "",
                      "END": "", "QSIZE": "",
                      "IDENTITY": "", "CHRO": "", "STRAND": "", "START.1": "", "END.1": "", "SPAN": "", "Exon": ""
                      }
        master_df = master_df.append(dummy_data, ignore_index=True)
        return master_df

    def _prepare_reduced_df(self):
        def create_pos19_str(x):
            return x["chromosome"] + ":" + str(x["start [hg19]"]) + "-" + str(x["end [hg19]"])

        self.reduced_df["position [hg19]"] = self.reduced_df.apply(create_pos19_str, axis=1)
        self.reduced_df["gene name"] = self.search_gene
        self.reduced_df["sequence length"] = self.reduced_df["end [hg19]"] - self.reduced_df["start [hg19]"]
        self.reduced_df["genomic sequence [hg19]"] = ""

    def _get_and_append_dna(self):
        pool = multiprocessing.Pool(processes=1)

        results_for_pseudogenes = []

        def apply_results(args):
            pseudogene_table, data, current_row = args
            # This is called whenever foo_pool(i) returns a result.
            # result_list is modified only by the main process, not the pool workers.
            self.reduced_df.set_value(current_row, "genomic sequence [hg19]", data)
            results_for_pseudogenes.append(pseudogene_table)

        result = pool.map(_get_data, ((i, row, self.hgsid, self.config, self.folder_name) for i, row in
                                      self.reduced_df.iterrows()))
        _ = [apply_results(x) for x in result if x is not None]
        pool.close()
        pool.join()

        return results_for_pseudogenes

    def _get_folder_name(self, gene_name):
        local_folder = os.path.join(os.getcwd(), gene_name)
        if os.path.isdir(local_folder):
            self._rename_old_data_folder(local_folder)

        os.makedirs(local_folder)
        return local_folder

    def _rename_old_data_folder(self, local_folder):
        msg = "Warning Folder {} already exists. The folder will be renamed. Continue? y or n".format(local_folder)
        if not self._unattended:
            if raw_input(msg).lower().strip() == "n":
                exit()
        else:
            back_fn_name = local_folder + "_" + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
            shutil.move(local_folder, back_fn_name)

    def find_pseudogenes(self):
        logger.info("Using hgsid: " + self.hgsid)
        logger.info("Source File: " + self.source_file)

        self.primers = self._prepare_input()

        if self.genes is None:
            print "Please enter a gene name"
            search_gene = raw_input('--> ')
            self.search_gene = search_gene.strip().upper()

            print "searching for: %s" % search_gene
            self.genes = search_gene
            if raw_input("Continue? y or n").lower().strip() == "n":
                exit()

        lines = self.genes.splitlines()
        for search_gene in lines:
            self.search_gene = search_gene
            self.reduced_df = self._find_genes()

            self.folder_name = self._get_folder_name(search_gene)

            self._prepare_reduced_df()
            self.results = self._get_and_append_dna()

            logger.info("Step 3: collecting results into output table")
            self._collect_results()

            logger.info("Step 4: generating DNA Table")
            file_name = "DNA FOR {}.xlsx".format(search_gene)

            local_result_path = os.path.join(self.folder_name, file_name)
            self.reduced_df.to_excel(local_result_path)
            logger.info(("Written to: {}".format(local_result_path)))
        logger.info(("Done"))


def _get_data(args):
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
            logger.info(data)
        except AttributeError:
            logger.info('No tables found, exiting')
            logger.warn(data)

        name = "{} {} {}-{}".format(row["gene name"], current_row, row["start [hg19]"], row["end [hg19]"])

        results_df = step2(data, name, cfg["detailed_results_to_excel"], folder)

        return results_df, data, current_row
    except Exception as ex:
        logger.exception("Something bad happened")
        logger.info(
            "Could not execute getDNA for {}, HGSID stale? Get a new one and paste it into config.json".format(
                dna_pos))


if __name__ == "__main__":
    worker = PseudoGeneFinder()
    worker.find_pseudogenes()
