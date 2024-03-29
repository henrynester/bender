import track, tline
import numpy as np, matplotlib.pyplot as plt
import ezdxf
from ezdxf import units
import math


class Bender:  # uh oh
    def __init__(self, track_sequence, floquet_unit_cell):
        self.track_sequence, self.floquet_unit_cell = track_sequence, floquet_unit_cell
        self.create_dxf()
        
    def construct_tline(self):
        self.unit_cell_vertices_count = self.floquet_unit_cell.vertices().shape[1]
        unit_cell_length = self.floquet_unit_cell.cell_length()
        arclength = trackseq.total_arclength()
        self.n_floquet_cells = int(arclength // unit_cell_length)-1
        
        tline_vertices = np.array([[], []])
        offset = np.array([[0.0], [0.0]])
        for i in range(self.n_floquet_cells):
            tline_vertices = np.append(tline_vertices, floquet.vertices() + offset, axis=1)
            offset += np.array([[unit_cell_length], [0.0]])

        self.tline_vertices = tline_vertices
        
    def bend_tline(self, plot=False):
        arclengths, heights = self.tline_vertices
        xs, ys, us, vs = self.track_sequence.vertices_normals(arclengths)

        self.bent_xs = xs + us * heights
        self.bent_ys = ys + vs * heights

        if plot:
            fig,ax=plt.subplots()
            ax.set_aspect('equal')
            ax.plot(self.bent_xs, self.bent_ys, marker='.', color='black')
            plt.show()

    def mirror_tline(self):
        self.tline_vertices[1,:] *= -1

    def create_dxf(self):
        self.doc = ezdxf.new()
        self.msp = self.doc.modelspace()
        self.doc.units = units.M

    def write_bent_tline(self):
        for n in range(self.n_floquet_cells):
            first=self.unit_cell_vertices_count*n
            last=first+self.unit_cell_vertices_count
            polygon = [(self.bent_xs[i], self.bent_ys[i]) for i in range(first,last)]
            self.msp.add_lwpolyline(polygon, close=True)

    def write_dxf(self,filename):
        # Save the DXF file
        self.doc.saveas(filename)

if __name__ == '__main__':
    fishboneA, fishboneB = tline.FishboneUnitCell(7e-6, 1e-6, 70e-6, 1e-6, gnd_spacing=2e-6, interdigitate=True), \
                            tline.FishboneUnitCell(7e-6, 1e-6, 100e-6, 1e-6, gnd_spacing=2e-6, interdigitate=True)
    floquet = tline.FloquetUnitCell()
    floquet.append_fishbones(fishboneA, 14)
    floquet.append_fishbones(fishboneB, 14)
    floquet.append_fishbones(fishboneA, 29)
    floquet.append_fishbones(fishboneB, 14)
    floquet.append_fishbones(fishboneA, 27)
    floquet.append_fishbones(fishboneB, 17)
    floquet.append_fishbones(fishboneA, 13)


    min_spacing = 4 * floquet.cell_min_spacing()
    turns = 2.9
    final_angle = turns * 2 * np.pi

    tline_compact_length = 25e-3

    fermat1 = track.FermatSpiralTrack(turns, min_spacing, True)
    fermat2 = track.FermatSpiralTrack(turns, min_spacing, False)
    arc2 = fermat2.construct_suitable_ArcTrack()
    arc1 = fermat1.construct_suitable_ArcTrack()
    straight2 = arc2.construct_suitable_StraightTrack(tline_compact_length)
    straight1 = arc1.construct_suitable_StraightTrack(tline_compact_length)

    trackseq = track.TrackSequence()
    trackseq.append_track(straight2)
    trackseq.append_track(arc2)
    trackseq.append_track(fermat2)
    trackseq.append_track(fermat1)
    trackseq.append_track(arc1)
    trackseq.append_track(straight1)

    bender = Bender(trackseq, floquet)
    bender.construct_tline()
    bender.bend_tline(plot=False)
    bender.write_bent_tline()

    bender.mirror_tline()
    bender.bend_tline(plot=True)
    bender.write_bent_tline()

    print(floquet.cell_length())
    print(trackseq.total_arclength())

    bender.write_dxf('output.dxf')
