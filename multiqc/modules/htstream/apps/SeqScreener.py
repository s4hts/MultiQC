from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import bargraph

#################################################

""" SeqScreener submodule for HTStream charts and graphs """

#################################################


class SeqScreener:

    ########################
    # Info about App
    def __init__(self):
        self.info = "A simple sequence screening tool which uses a kmer lookup approach to identify reads from an unwanted source."
        self.type = "read_reducer"


    ########################
    # Bargraph Function
    def bargraph(self, json, reads_screened, index):

        # configuration dictionary for bar graph
        config = {
            "title": "HTStream: Sequences Identified",
            "id": "htstream_seqscreener_bargraph_" + index,
            "ylab": "Reads",
        }

        html = ""

        # if no overlaps at all are present, return nothing
        if reads_screened == 0:
            html += '<div class="alert alert-info"> <strong>Notice:</strong> No reads were identified from samples. </div>'
            return html

        # bargraph dictionary. Exact use of example in MultiQC docs.
        categories = OrderedDict()

        # Colors for sections
        categories["Ss_PE_hits" + index] = {"name": "Paired End", "color": "#779BCC"}
        categories["Ss_SE_hits" + index] = {"name": "Single End", "color": "#D1ADC3"}

        # create bargrpah
        html += bargraph.plot(json, categories, config)

        return html

    ########################
    # Main Function
    def execute(self, json, index):

        stats_json = OrderedDict()
        overview_dict = {}
        reads_screened = 0 

        for key in json.keys():

            pe_hits = json[key]["Paired_end"]["hits"]    
            se_hits = json[key]["Single_end"]["hits"]

            reads_screened += pe_hits + se_hits

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Reads_Lost": (json[key]["Fragment"]["in"] - json[key]["Fragment"]["out"])
                / json[key]["Fragment"]["in"],
                "Percent_Hits": (pe_hits + se_hits) / json[key]["Fragment"]["in"],
            }

            # sample entry for stats dictionary
            stats_json[key] = {
                "Ss_PE_hits" + index: pe_hits,
                "Ss_SE_hits" + index: se_hits,
            }

        # sections and figure function calls
        section = {"Bargraph": self.bargraph(stats_json, reads_screened, index), "Overview": overview_dict}

        return section
