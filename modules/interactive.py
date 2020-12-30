#!/usr/bin/env python3

from . import common, annotate
import statistics
from Bio import SeqIO
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import plot

# this section silences warnings while importing an improperly formatted genbank file.
# Shouldn't need when program is called to plot from the annotate module since it will pass a SeqRecord
# instead of a text file
import warnings
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)


class Entry:
    """Entry converts a SeqFeature into an object that has all info needed for plotting a single data point."""
    def __init__(self, args, seqfeature, normalized=False):
        self.args = args
        self.feature = seqfeature
        self.contig = self.feature.qualifiers['contig_id'][0]
        self.feature_type = self.retrieve_feature_type()
        self.normalized = normalized
        self.position = self.feature.location
        self.start = int(self.feature.location.start)
        self.annotation = self.retrieve_gene_annotation()
        self.locus = self.feature.qualifiers['locus_tag'][0]
        self.length = len(self.feature)
        self.blast_lengths = self.retrieve_blasthit_lengths()
        self.mean_blast_length = self.retrieve_mean_blast_length()
        self.blast_stdev = self.blast_stdev()

        if normalized:
            self.length = self.normalize_length()
            self.blast_stdev = self.normalize_stdev()   # must be called before normalizing length
            self.mean_blast_length = 1
            self.normalized_length = self.normalize_length()

    def __str__(self):
        return ", ".join([self.locus, self.feature_type, str(self.position), str(self.length),
                          "%s blast hits" % len(self.blast_lengths)])

    def retrieve_blasthit_lengths(self):
        lengths = []
        try:
            for hit in self.feature.qualifiers['hits']:
                if isinstance(hit, str):
                    hit = convert_string_to_blasthit(hit)
                lengths.append(annotate.blasthit_length(hit))
        except KeyError:
            pass
        return lengths

    def retrieve_feature_type(self):
        if self.feature.type == 'pseudogene':
            if isinstance(self.feature.qualifiers['pseudo_type'], annotate.PseudoType):
                return str(self.feature.qualifiers['pseudo_type'])
            else:
                return str(self.feature.qualifiers['pseudo_type'][0])
        else:
            return self.feature.type

    def retrieve_gene_annotation(self):
        if self.feature_type == 'intergenic':
            annotation = 'Intergenic'
        elif self.feature_type == 'PseudoType.intergenic':
            annotation = 'Intergenic pseudo'
        else:
            try:
                product = self.feature.qualifiers['product'][0]
            except KeyError:
                product = 'unknown product'
            try:
                gene = self.feature.qualifiers['gene'][0]
            except KeyError:
                gene = 'unknown gene'
            annotation = "%s, %s" % (gene, product)
        return annotation

    def retrieve_mean_blast_length(self):
        if len(self.blast_lengths) == 0:
            return 0
        else:
            return sum(self.blast_lengths) / len(self.blast_lengths)

    def blast_stdev(self):
        if len(self.blast_lengths) < 2:
            return 0
        else:
            return statistics.stdev(self.blast_lengths)

    def normalize_length(self):
        if len(self.blast_lengths) == 0:
            return 1
        else:
            return self.length / self.mean_blast_length

    def normalize_stdev(self):
        if len(self.blast_lengths) < 2:
            return 0
        else:
            return statistics.stdev(self.blast_lengths) / self.mean_blast_length

    def normalize(self):
        """Returns a new copy of itself with normalized values"""
        return Entry(args=self.args,
                     seqfeature=self.feature,
                     normalized=True)


class Figure:
    """Calling this object will create a figure with the provided datasets."""
    def __init__(self, args, file_dict, datasets, descriptions):
        self.args = args
        self.file_dict = file_dict
        self.datasets = datasets
        self.num_plots = len(datasets)
        self.fig = make_subplots(rows=self.num_plots, cols=1, subplot_titles=descriptions)

        for i, dataset in enumerate(self.datasets):
            self.generate_subplot(dataset, i+1)
        self.fig.update_layout(hovermode='x unified',
                               xaxis=dict(ticks="",
                                          showticklabels=False))

    def show(self):
        self.fig.show()

    def offline_plot(self):
        plot(self.fig, filename=self.file_dict['interactive'], auto_open=False)

    def bar_colours(self, dataset):
        colours = []
        for entry in dataset:
            if "Pseudo" in str(entry.feature_type):
                colours.append("#E5446D")
            elif entry.feature_type == 'intergenic':
                colours.append("#BAC2C9")
            elif len(entry.blast_lengths) == 0:
                colours.append("#465362")
            else:
                colours.append("#70A0FF")
        return colours

    def hovertext(self, dataset, object_type):
        hovertext = []
        for entry in dataset:
            if object_type == "bar":
                bar_text = "Contig: %s %s<br>" \
                           "Locus: %s<br>" \
                           "Annotation: %s<br>" % (entry.contig, str(entry.position), entry.locus, entry.annotation)
                if 'pseudo' in entry.feature_type.lower():
                    bar_text = bar_text + "Pseudo call: %s" % entry.feature_type
                hovertext.append(bar_text)
            elif object_type == "scatter":
                hovertext.append("%s blast hits." % len(entry.blast_lengths))
            else:
                raise RuntimeError("wrong use of Figure.hovertext().")
        return hovertext

    def trace_name(self, dataset, object_type):
        if object_type == 'bar':
            if data_is_normalized(dataset):
                trace_name = 'Normalized length'
            else:
                trace_name = 'Length (nt)'

        elif object_type == 'scatter':
            if data_is_normalized(dataset):
                trace_name = 'Normalized BLAST hit length'
            else:
                trace_name = 'BLAST mean hit length (nt)'

        else:
            raise RuntimeError("wrong use of Figure.trace_name().")
        return trace_name

    def generate_subplot(self, dataset, row: int):
        x_vals = np.arange(len(dataset))
        y_mean_hit_vals = [x.mean_blast_length for x in dataset]
        y_mean_hit_stdev = [x.blast_stdev for x in dataset]
        y_gene_vals = [x.length for x in dataset]
        colours = self.bar_colours(dataset)
        markersize = 1
        error = dict(type='data',
                     array=y_mean_hit_stdev,
                     width=markersize,
                     thickness=markersize,
                     visible=True)

        self.fig.add_trace(go.Bar(x=x_vals,
                                  y=y_gene_vals,
                                  marker=dict(color=colours),
                                  hoverlabel=dict(namelength=-1),
                                  hovertext=self.hovertext(dataset, 'bar'),
                                  showlegend=False,
                                  name=self.trace_name(dataset, 'bar')),
                           row=row,
                           col=1)

        if data_is_normalized(dataset):
            self.fig.add_trace(go.Scatter(x=x_vals,
                                          y=y_mean_hit_vals,
                                          error_y=error,
                                          mode='lines',
                                          marker=dict(color='black', size=0),
                                          line=dict(width=markersize),
                                          hoverlabel=dict(namelength=-1),
                                          hovertext=self.hovertext(dataset, 'scatter'),
                                          showlegend=False,
                                          name=self.trace_name(dataset, 'scatter')),
                               row=row,
                               col=1)
        else:
            self.fig.add_trace(go.Scatter(x=x_vals,
                                          y=y_mean_hit_vals,
                                          error_y=error,
                                          mode='markers',
                                          marker=dict(color='black', size=markersize),
                                          hoverlabel=dict(namelength=-1),
                                          hovertext=self.hovertext(dataset, 'scatter'),
                                          showlegend=False,
                                          name=self.trace_name(dataset, 'scatter')),
                               row=row,
                               col=1)

        self.fig.update_xaxes(ticks="", showticklabels=False, row=row, col=1)


def data_is_normalized(dataset):
    """
    Returns true if every entry in the dataset is normalized.
    """
    return all([entry.normalized for entry in dataset])


def gbk_to_seqrecord_list(args, gbk: str):
    """
    Uses biopython SeqIO to convert a genbank file into a list of SeqRecord, one for each contig
    """
    seqrecord_list = []
    with open(gbk, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            seqrecord_list.append(record)
    return seqrecord_list


def convert_string_to_blasthit(string):
    fields = string.replace("BlastHit(", "")[:-1].split(",")
    fields_dict = {}
    for field in fields:
        if 'stitle' not in field:
            try:
                key, value = [x.strip(" ") for x in field.split("=")]
                fields_dict[key] = common.literal_eval(value.strip("\'"))
            except ValueError:
                pass
        fields_dict['stitle'] = "info lost"
    blasthit = annotate.BlastHit(**fields_dict)
    return blasthit


def retrieve_blasthit_lengths(args, hits):
    lengths = []
    for hit in hits:
        if isinstance(hit, str):
            hit = convert_string_to_blasthit(hit)
        lengths.append(annotate.blasthit_length(hit))
    return lengths


def features_to_data(args, features):
    data = []
    for feature in features:
        data.append(Entry(args=args, seqfeature=feature))
    return data


def genome_to_graphs(args, file_dict, genome, show=False):
    # Collect all relevant features
    cds_features = annotate.extract_features_from_genome(args, genome, 'CDS')
    intergenic_features = annotate.extract_features_from_genome(args, genome, 'intergenic')
    pseudo_features = annotate.extract_features_from_genome(args, genome, 'pseudogene')
    all_features = cds_features + intergenic_features + pseudo_features

    # Convert to data (list of Entry objects)
    data = features_to_data(args, all_features)
    data.sort(key=lambda x: x.start)
    data.sort(key=lambda x: x.contig)  # Result is data sorted by contig, then position on each contig

    # Each description and data will generate one plot
    desc1 = "All regions sorted by contig then location on contig"
    plot1_data = data

    desc2 = "Intergenic regions removed, gene lengths normalized with mean blast hit lengths"
    plot2_data = [d.normalize() for d in data if d.feature_type != 'intergenic']

    desc3 = "Sorted by normalized gene length"
    plot3_data = sorted(plot2_data, key=lambda x: x.length)

    # Plot the data
    fig = Figure(args=args,
                 file_dict=file_dict,
                 datasets=[plot1_data, plot2_data, plot3_data],
                 descriptions=[desc1, desc2, desc3])

    if show:
        fig.show()
    else:
        fig.offline_plot()


def main():
    args = common.get_args('interactive')
    annotated_genome = gbk_to_seqrecord_list(args, args.annotated_genome)
    if args.outprefix:
        show = False
        file_dict = common.file_dict(args)
    else:
        show = True
        file_dict = common.file_dict(args, outprefix='Irrelevant')
    genome_to_graphs(args, file_dict, annotated_genome, show)


if __name__ == '__main__':
    main()
