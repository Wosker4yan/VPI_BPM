"""
Created on T Sep 30 12:01 2022

Class to get the polygons from the VPI
modified slightly by Vahram

Reusing for VPI


@author: Marco
BRIGHT PHOTONICS
"""

import nazca as nd

class Masktopolygon:
    def __init__(self, layer_dic=None, material_index=None):
        if layer_dic is None:
            xs = nd.add_xsection('VPI')
            nd.add_layer2xsection('VPI', 1)

            layer_dic = {
                1: 2.0,
            }
        self.layer_dic = {nd.get_layer(lay): val for lay, val in layer_dic.items()}

    def get_default_ic(self, width=2, radius=20.0):
        xs = nd.get_xsection('VPI')
        return nd.interconnects.Interconnect(xs='VPI', width=width, radius=radius)

    def get_polygons(self, cell):
        polygon_dic = {lay: [] for lay in self.layer_dic}

        for param in nd.cell_iter(cell, revisit=True):
            if param.cell_open:
                pointer = param.transflip_glob[0]
                for poly, xy, bb in param.iters['polygon']:
                    if poly.layer in self.layer_dic:
                        points = [pointer.copy().move(*x).xy() for x in xy]
                        polygon_dic[poly.layer].append(points)
        return polygon_dic

