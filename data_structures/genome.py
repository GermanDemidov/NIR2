# -*- coding: utf-8 -*-
from data_structures.chromosome import Chromosome

__author__ = "Sergey Aganezov"
__email__ = "aganezov@gwu.edu"
__status__ = "develop"


class Genome(object):
    def __init__(self, name, manager, resolution):
        self.name = name
        self.chromosomes = {}
        self.location_info = {}
        self.mapping = {}
        self.manager = manager
        self.resolution = resolution

    def get_chromosome_by_sb(self, sb):
        if sb in self.location_info:
            return self.location_info[sb]

    def add_sb(self, chromosome_name, sb):
        if chromosome_name not in self.chromosomes:
            self._create_and_add_chromosome(chromosome_name)
        chromosome = self.chromosomes[chromosome_name]
        chromosome.add_sb(sb)

    def _add_chromosome(self, chromosome):
        self.chromosomes[chromosome] = chromosome

    def _create_and_add_chromosome(self, chromosome_name):
        chromosome = Chromosome(name=chromosome_name, genome=self)
        self._add_chromosome(chromosome)

    def add_mapping(self, res, sb_name, mapped_sequence):
        chromosome = self.location_info[sb_name]
        chromosome.add_mapping(sb_name, res, mapped_sequence)

    def retrieve_sb_by_name(self, *args):
        result = []
        for entry in args:
            chromosome = self.location_info[entry]
            result.append(chromosome.sb_name_mapping[entry])
        return result

    def get_mapping_and_offset_for_sb(self, sb, res_to_map):
        chromosome = self.location_info[sb]
        sb = chromosome.sb_name_mapping[sb]
        l_map, r_map, l_offset, r_offset = chromosome.get_two_closest_mappings_and_relative_position(sb=sb,
                                                                                                     res_to_map=res_to_map)
        if l_map == r_map:
            mapping_start = chromosome.get_left_outer_most_mapping(sb=l_map, res_to_map=res_to_map)
        else:
            mapping_start = chromosome.get_right_outer_most_mapping(sb=l_map, res_to_map=res_to_map)
        mapping_end = chromosome.get_left_outer_most_mapping(sb=r_map, res_to_map=res_to_map)
        return mapping_start, mapping_end, l_offset, r_offset

    def get_mapped_set(self, sb_start, sb_end, l_offset=0, r_offset=None):
        chromosome1, chromosome2 = self.location_info[sb_start], self.location_info[sb_end]
        assert chromosome1 == chromosome2
        vague_mapping = chromosome1.get_sbs_between(sb_start=sb_start, sb_end=sb_end)
        return set(vague_mapping[l_offset:r_offset])