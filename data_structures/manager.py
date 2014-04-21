# -*- coding: utf-8 -*-
from collections import defaultdict
from data_structures.genome import Genome

__author__ = 'Sergey Aganezov'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'


class Manager(object):
    def __init__(self):
        self.resolutions = defaultdict(dict)

    def _add_genome(self, genome, resolution_name):
        self.resolutions[resolution_name][genome] = genome

    def _create_and_add_genome(self, genome_name, resolution_name):
        genome = Genome(name=genome_name, manager=self, resolution=resolution_name)
        self._add_genome(genome=genome, resolution_name=resolution_name)

    def add_sb(self, resolution_name, genome_name, chromosome_name, sb_info):
        if genome_name not in self.resolutions[resolution_name]:
            self._create_and_add_genome(genome_name=genome_name, resolution_name=resolution_name)
        genome = self.resolutions[resolution_name][genome_name]
        genome.add_sb(chromosome_name, sb_info)

    def get_mapping_set_for_sb(self, genome_name, sb_name, res1, res2):
        if genome_name in self.resolutions[res1] and genome_name in self.resolutions[res2]:
            genome = self.resolutions[res1][genome_name]
            mapped_genome = self.resolutions[res2][genome_name]
            l_map, r_map, l_offset, r_offset = genome.get_mapping_and_offset_for_sb(sb_name, res2)
            if l_map == r_map:
                return set()
            return mapped_genome.get_mapped_set(sb_start=l_map, sb_end=r_map, l_offset=l_offset, r_offset=r_offset)



