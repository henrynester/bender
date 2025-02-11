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
    filename = 'test1.dxf'
    feature_size = 1.5e-6
    wcenter = feature_size
    wline = 45e-6 - wcenter
    wload = 30e-6 - wcenter
    fishboneA, fishboneB = tline.FishboneUnitCell(8e-6, feature_size, wline/2, wcenter, gnd_spacing=2e-6, interdigitate=True), \
                            tline.FishboneUnitCell(8e-6, feature_size, wload/2, wcenter, gnd_spacing=2e-6, interdigitate=True)
    # wcenter = 8e-6
    # wline = 40e-6 - wcenter
    # fishboneA = tline.FishboneUnitCell(4e-6, 2e-6, wline/2, wcenter, gnd_spacing=2e-6, interdigitate=False)
    floquet = tline.FloquetUnitCell()
    floquet.append_fishbones(fishboneA, 1)





    #Alec 10GHz
    # floquet.append_fishbones(fishboneA, 19)
    # floquet.append_fishbones(fishboneB, 8)
    # floquet.append_fishbones(fishboneA, 38)
    # floquet.append_fishbones(fishboneB, 8)
    # floquet.append_fishbones(fishboneA, 37)
    # floquet.append_fishbones(fishboneB, 10)
    # floquet.append_fishbones(fishboneA, 18)

    #Jordan 5GHz
    # floquet.append_fishbones(fishboneA, 39)
    # floquet.append_fishbones(fishboneB, 13)
    # floquet.append_fishbones(fishboneA, 79)
    # floquet.append_fishbones(fishboneB, 13)
    # floquet.append_fishbones(fishboneA, 78)
    # floquet.append_fishbones(fishboneB, 14)
    # floquet.append_fishbones(fishboneA, 39)


    min_spacing = 5.5*floquet.cell_min_spacing()
    turns = 4
    final_angle = turns * 2 * np.pi

    tline_compact_length = 100e-5

    origin = np.array((0, 0))
    end = np.array((tline_compact_length, 0))
    straight1 = track.StraightTrack(origin, end)
    r = 20e-5
    arc1 = track.ArcTrack(np.array((tline_compact_length, r)),r,-np.pi/2,np.pi/2, False)
    straight2 = track.StraightTrack(np.array((tline_compact_length, 2*r)), np.array((0, 2*r)))

    # fermat1 = track.FermatSpiralTrack(turns, min_spacing, True)
    # fermat2 = track.FermatSpiralTrack(turns, min_spacing, False)
    # arc2 = fermat2.construct_suitable_ArcTrack()
    # arc1 = fermat1.construct_suitable_ArcTrack()
    # straight2 = arc2.construct_suitable_StraightTrack(tline_compact_length)
    # straight1 = arc1.construct_suitable_StraightTrack(tline_compact_length)

    trackseq = track.TrackSequence()
    # trackseq.append_track(straight2)
    # trackseq.append_track(arc2)
    # trackseq.append_track(fermat2)
    # trackseq.append_track(fermat1)
    trackseq.append_track(straight1)
    trackseq.append_track(arc1)
    trackseq.append_track(straight2)

    bender = Bender(trackseq, floquet)
    bender.construct_tline()
    bender.bend_tline(plot=True)
    bender.write_bent_tline()

    # bender.mirror_tline()
    # bender.bend_tline(plot=True)
    # bender.write_bent_tline()

    print(floquet.cell_length())
    print(trackseq.total_arclength())
    print(bender.n_floquet_cells)

    bender.write_dxf(filename)

