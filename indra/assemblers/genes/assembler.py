# -*- coding: utf-8 -*-

"""An assembler for gene lists and related gene-centric outputs."""

from __future__ import absolute_import, print_function, unicode_literals

import csv
import os

from jinja2 import Template

__all__ = [
    'GeneAssembler',
]

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)))

#: Path to RefSeq data for Ideogram output
REFSEQ_PATH = os.path.join(HERE, 'refseq_human.csv')

TEMPLATE_PATH = os.path.join(HERE, 'ideogram_template.html')
with open(TEMPLATE_PATH, 'rt') as f:
    html_template = Template(f.read())


class GeneAssembler:
    """The Gene assembler assembles INDRA Statements into a gene list.

    This graph can then be used with the GSEA software from the Broad
    Institute, or output as a visualization with Ideogram.

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be added to the assembler's list
        of Statements.

    Attributes
    ----------
    genes : set[str]
        A set of strings representing gene names in the model.
    """

    def __init__(self, stmts=None):
        self.stmts = stmts or []
        self.genes = set()

    def make_model(self):
        """Assemble the graph from the assembler's list of INDRA Statements."""
        for statement in self.stmts:
            try:
                agents = statement.agent_list()
            except Exception:
                continue
            else:
                for agent in agents:
                    if agent is None:
                        continue
                    if 'HGNC' in agent.db_refs:
                        self.genes.add(agent.name)

        return self.genes

    def save_model(self, path, format='gsea'):
        """Save the assembled model's SIF string into a file.

        Parameters
        ----------
        path : str
            The name of the file to save the SIF into.
        """
        with open(path, 'w') as file:
            if format == 'gsea':
                print(self.to_gsea_str(), file=file)
            elif format == 'ideogram':
                print(self.to_ideogram_html(), file=file)

    def to_gsea_str(self, first='# INDRA Genes'):
        """Return a SIF string of the assembled model."""
        return first + '\n' + '\n'.join(sorted(self.genes))

    def to_ideogram_html(self):
        """Get an Ideogram HTML document as a string."""
        return html_template.render(
            annotations=self.get_ideogram_annotations(),
        )

    def get_ideogram_annotations(self):
        with open(REFSEQ_PATH) as file:
            with csv.reader(file) as reader:
                next(reader)  # skip header
                return [
                    {
                        'name': symbol,
                        'start': int(start),
                        'stop': int(stop),
                    }
                    for _, symbol, start, stop in reader
                    if symbol not in self.genes
                ]
