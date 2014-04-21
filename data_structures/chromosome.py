# -*- coding: utf-8 -*-
from collections import defaultdict
from data_structures.sb import StartTelomereSyntenyBlock, EndTelomereSyntenyBlock, SyntenyBlock
from bintrees import AVLTree

__author__ = 'Sergey Aganezov'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'


class Chromosome(object):
    def __init__(self, name, genome):
        self.name = name
        self.genome = genome
        self.sb_name_mapping = {}
        self.start_telomer = StartTelomereSyntenyBlock()
        self.end_telomer = EndTelomereSyntenyBlock()
        self.blocks = AVLTree()
        self.add_sb(self.start_telomer)
        self.add_sb(self.end_telomer)
        # maps sb into a dict of resolutions, it has mapping in
        # in each particular resolution, mapping is represented as a list,
        # this particular block is mapped to
        self.mapping = defaultdict(lambda: defaultdict(list))

    def add_sb(self, sb):
        if not type(sb, SyntenyBlock):
            sb = self._create_sb(sb)
        self.genome.location_info[sb] = self
        self.sb_name_mapping[sb] = sb
        self.blocks.insert(key=(sb.start, sb.end), value=sb)

    def _create_sb(self, info):
        sb = SyntenyBlock(name=info.name,
                          start=info.start,
                          end=info.end,
                          strand=info.strand,
                          chromosome=self)
        return sb

    def add_mapping(self, sb, res_to_map, sbs):
        for synteny_block in sbs:
            if synteny_block not in self.mapping[sb][res_to_map]:
                self.mapping[sb][res_to_map].append(synteny_block)

    def get_right_outer_most_mapping(self, sb, mapping):
        if sb in self.mapping and mapping in self.mapping[sb]:
            return self.mapping[sb][mapping][-1]
        return self.start_telomer

    def get_left_outer_most_mapping(self, sb, mapping):
        if sb in self.mapping and mapping in self.mapping[sb]:
            return self.mapping[sb][mapping][0]
        return self.end_telomer

    def get_sbs_between(self, sb_start, sb_end):
        # slice from 1: because, value_slice returns all values, which keys are greater or equals that minimum,
        # but less than maximum. Because of that "or equal" we have to slice
        result = list(self.blocks.value_slice((sb_start.start, sb_start.end), (sb_end.start, sb_end.end)))
        return result[1:] if len(result) > 0 else result

    def get_two_closest_mappings_and_relative_position(self, sb, res_to_map):
        l_position, r_position = 0, 0
        if res_to_map in self.mapping[sb]:
            return sb, sb, l_position, r_position
        for value in self.blocks.value_slice(start_key=(sb.start, sb.end)):
            r_position += 1
            if res_to_map in self.mapping[value]:
                right_closest_mapping = value
                break
        else:
            right_closest_mapping = self.end_telomer
        for value in self.blocks.value_slice(end_key=(sb.start, sb.end),
                                             reverse=True):
            l_position += 1
            if res_to_map in self.mapping[value]:
                left_closest_mapping = value
                break
        else:
            left_closest_mapping = self.start_telomer
        return left_closest_mapping, right_closest_mapping, l_position, r_position

    def __hash__(self):
        return self.name.__hash__()

