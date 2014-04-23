# -*- coding: utf-8 -*-

__author__ = 'Sergey Aganezov'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'


class SyntenyBlock(object):
    def __init__(self, name, start, end, strand, chromosome=None):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.chromosome = chromosome

    def __hash__(self):
        return self.name.__hash__()

    @property
    def genome(self):
        return self.chromosome.genome


class StartTelomereSyntenyBlock(SyntenyBlock):
    def __init__(self, chromosome_name):
        super(StartTelomereSyntenyBlock, self).__init__(name="start_telomere",
                                                        start=-1, end=0,
                                                        strand="+")
        self.name += chromosome_name


class EndTelomereSyntenyBlock(SyntenyBlock):
    def __init__(self, chromosome_name):
        super(EndTelomereSyntenyBlock, self).__init__(name="end_telomere",
                                                      start=float("inf"),
                                                      end=float("inf"),
                                                      strand="+")
        self.name += chromosome_name
