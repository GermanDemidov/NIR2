# -*- coding: utf-8 -*-
__author__ = "Sergey Aganezov"
__email__ = "aganezov@gwu.edu"
__status__ = "develop"


class Genome(object):
    def __init__(self):
        self.chromosomes = {}
        self.location_info = {}
        self.mapping = {}

    def get_chromosome_by_sb(self, sb):
        if sb in self.location_info:
            return self.location_info[sb]

    def add_sb(self, chromosome, sb):
        if chromosome not in self.chromosomes:
            self.add_chromosome(chromosome)
        self.location_info[sb] = chromosome
        chromosome.add_sb(sb)

    def add_chromosome(self, chromosome):
        self.chromosomes[chromosome] = chromosome