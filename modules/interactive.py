#!/usr/bin/env python3

from . import common, annotate
import statistics
from Bio import SeqIO
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import plot


class Entry:
    """Entry converts a SeqFeature into an object that has all info needed for plotting a single data point."""
    def __init__(self, args, seqfeature, normalized=False):
        self.args = args
        self.feature = seqfeature
        if isinstance(self.feature.qualifiers['contig_id'], list):
            self.contig = self.feature.qualifiers['contig_id'][0]
        else:
            self.contig = self.feature.qualifiers['contig_id']
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
    def __init__(self, args, file_dict):
        self.args = args,
        self.file_dict = file_dict
        self.fig = None
        self.config = self.plot_config()
        self.filename = None

    def plot_config(self):
        config = dict(displayModeBar=True,
                      displaylogo=False)
        return config

    def show(self):
        self.fig.show(config=self.config)

    def offline_plot(self):
        plot(self.fig, filename=self.filename, auto_open=False)

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
            bar_text = "Contig: %s<br>" \
                       "Locus: %s %s<br>" \
                       "Annotation: %s<br>" % (entry.contig, entry.locus, str(entry.position), entry.annotation)

            if 'pseudo' in entry.feature_type.lower():
                bar_text = bar_text + "Pseudo call: %s" % entry.feature_type

            scatter_text = "%s blast hits." % len(entry.blast_lengths)
            pie_length = 'Length (nt): %s<br>' \
                         'BLAST hit length (Mean nt, SD): %s +/- %s<br>' % (entry.length,
                                                                            round(entry.mean_blast_length, 3),
                                                                            round(entry.blast_stdev, 3))
            pie_text = pie_length + bar_text + "<br>" + scatter_text

            if object_type == "bar":
                hovertext.append(bar_text)
            elif object_type == "scatter":
                hovertext.append(scatter_text)
            elif object_type == 'pie':
                hovertext.append(pie_text)
            else:
                raise RuntimeError("wrong use of Figure.hovertext().")

        return hovertext

    def show_alignment(self):
        """This function should bring up the associated protein alignment."""
        return False


class Bar(Figure):
    """Calling this object will create a figure with the provided datasets."""
    def __init__(self, args, file_dict, datasets, descriptions):
        super().__init__(args, file_dict)
        self.filename = self.file_dict['interactive_bar']
        self.datasets = datasets
        self.num_plots = len(datasets)
        self.fig = make_subplots(rows=self.num_plots, cols=1, subplot_titles=descriptions)
        for i, dataset in enumerate(self.datasets):
            self.generate_subplot(dataset, i+1)

        self.fig.update_layout(hovermode='x unified',
                               xaxis=dict(ticks="",
                                          showticklabels=False),
                               width=self.fig_width())

    def fig_width(self):
        """In pixels. Below 350 features, lets plotly set this automatically."""
        largest_dataset = sorted(self.datasets, key=lambda x: len(x))[-1]
        dataset_length = len(largest_dataset)
        if dataset_length < 350:
            return None
        else:
            return dataset_length * 5

    def trace_name(self, dataset, object_type):
        if object_type == 'bar':
            if data_is_normalized(dataset):
                trace_name = 'Normalized length'
            else:
                trace_name = 'Length (nt)'

        elif object_type == 'scatter':
            if data_is_normalized(dataset):
                trace_name = 'Normalized BLAST hit length (Mean, SD)'
            else:
                trace_name = 'BLAST hit length (Mean nt, SD)'

        else:
            raise RuntimeError("wrong use of Figure.trace_name().")
        return trace_name

    def generate_subplot(self, dataset, row: int):
        # Collect all the data from the dataset
        x_vals = np.arange(len(dataset))
        y_mean_hit_vals = [x.mean_blast_length for x in dataset]
        y_mean_hit_stdev = [x.blast_stdev for x in dataset]
        y_gene_vals = [x.length for x in dataset]

        # Any common keyword arguments for every graph_object go here
        markersize = 1
        go_kwargs = dict(x=x_vals,
                         hoverlabel=dict(namelength=-1),
                         showlegend=False)
        # Error bars
        error = dict(type='data',
                     array=y_mean_hit_stdev,
                     width=markersize,
                     thickness=markersize,
                     visible=True)
        # Bar plot
        bar = go.Bar(go_kwargs,
                     y=y_gene_vals,
                     marker=dict(color=self.bar_colours(dataset)),
                     hovertext=self.hovertext(dataset, 'bar'),
                     name=self.trace_name(dataset, 'bar'))

        # If normalized to y=1, no need for a marker but a line will be drawn at y=1.
        if data_is_normalized(dataset):
            scatter = go.Scatter(go_kwargs,
                                 y=y_mean_hit_vals,
                                 error_y=error,
                                 mode='lines',
                                 marker=dict(color='black', size=0),
                                 line=dict(width=markersize),
                                 hovertext=self.hovertext(dataset, 'scatter'),
                                 name=self.trace_name(dataset, 'scatter'))
        # If not normalized, no line but draw markers for each point
        else:
            scatter = go.Scatter(go_kwargs,
                                 y=y_mean_hit_vals,
                                 error_y=error,
                                 mode='markers',
                                 marker=dict(color='black', size=markersize),
                                 hovertext=self.hovertext(dataset, 'scatter'),
                                 name=self.trace_name(dataset, 'scatter'))

        self.fig.add_trace(bar, row=row, col=1)
        self.fig.add_trace(scatter, row=row, col=1)
        self.fig.update_xaxes(range=[-1, x_vals[-1] + 1], ticks="", showticklabels=False, row=row, col=1)
        self.fig.update_yaxes(fixedrange=True, row=row, col=1)
        self.fig.layout.annotations[row-1].update(xanchor='left', x=0)


class Map(Figure):
    """Makes a genome diagram from the plotly piechart"""
    def __init__(self, args, file_dict, dataset):
        super().__init__(args, file_dict)
        self.filename = self.file_dict['interactive_map']
        self.dataset = dataset
        self.data_by_contig = self.split_dataset_by_contig()
        self.num_plots = len(self.data_by_contig)
        self.fig = make_subplots(cols=self.num_plots, rows=1, specs=[[{'type': 'domain'}]*self.num_plots])
        for i, dataset in enumerate(self.data_by_contig):
            self.genome_map(dataset['contig'], dataset['entries'], i + 1)

    def split_dataset_by_contig(self):
        replicon_limit = 5
        datasets = []
        contigs = sorted(list(set([x.contig for x in self.dataset])))

        for contig in contigs:
            datasets.append({'contig': contig,
                             'entries': self.entries_on_contig(contig)})

        if len(datasets) >= replicon_limit:
            common.print_with_time('%s replicons in input file. Will only display visualization of first %s replicons.' % (len(datasets), replicon_limit))
            datasets = datasets[:replicon_limit]

        return datasets

    def entries_on_contig(self, contig):
        entries = []

        for entry in self.dataset:
            if entry.contig == contig:
                entries.append(entry)

        return entries

    def genome_map(self, contig, dataset, col: int):
        """
        Makes a genome diagram using the plotly.pie object.
        TODO:
        - fix scaling for larger genomes so that they are readable
        """
        values = [entry.length for entry in dataset]
        labels = [entry.locus for entry in dataset]
        pie_kwargs = dict(direction='clockwise',
                          hole=0.90,
                          sort=False,
                          hoverinfo='text',
                          insidetextorientation='radial',
                          textinfo='none',
                          showlegend=False)
        pie = go.Pie(pie_kwargs,
                     values=values,
                     labels=labels,
                     title=contig,
                     marker=dict(colors=self.bar_colours(dataset),
                                 line=dict(color='black', width=0.01)),
                     hovertext=self.hovertext(dataset, 'pie'),
                     scalegroup='one')

        self.fig.add_trace(pie, row=1, col=col)


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
    all_features = annotate.extract_features_from_genome(args, genome, ['CDS', 'intergenic', 'pseudogene'])

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
    bar = Bar(args=args,
              file_dict=file_dict,
              datasets=[plot1_data, plot2_data, plot3_data],
              descriptions=[desc1, desc2, desc3])

    map = Map(args=args,
              file_dict=file_dict,
              dataset=plot1_data)
    if show:
         bar.show()
         #map.show()
    else:
        bar.offline_plot()
        #map.offline_plot()


def main():
    # this import silences warnings while importing an improperly formatted genbank file.
    # Shouldn't need when program is called to plot from the annotate module since it will pass a SeqRecord
    # instead of a text file
    import warnings
    from Bio import BiopythonParserWarning
    warnings.simplefilter('ignore', BiopythonParserWarning)

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
