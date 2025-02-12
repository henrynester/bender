import track, tline
import numpy as np, matplotlib.pyplot as plt
import ezdxf
from ezdxf import units
import ezdxf.path
import progress.bar
import json, sys

class Bender:  # uh oh   
    def __init__(self, track_sequence, floquet_unit_cell, fillet_radius=None):
        self.track_sequence, self.floquet_unit_cell = track_sequence, floquet_unit_cell
        self.fillet_radius = fillet_radius

        self.create_dxf()
        
    def construct_tline(self):
        self.unit_cell_vertices_count = self.floquet_unit_cell.vertices().shape[1]
        unit_cell_length = self.floquet_unit_cell.cell_length()
        arclength = trackseq.total_arclength()
        self.n_floquet_cells = int(arclength // unit_cell_length)

        bar = progress.bar.Bar('', max=self.n_floquet_cells)
        
        tline_vertices = np.array([[], []])
        offset = np.array([[0.0], [0.0]])
        for i in range(self.n_floquet_cells):
            tline_vertices = np.append(tline_vertices, floquet.vertices() + offset, axis=1)
            offset += np.array([[unit_cell_length], [0.0]])
            bar.next()
        bar.finish()

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
        bar = progress.bar.Bar('', max=self.n_floquet_cells)
        for n in range(self.n_floquet_cells):
            first=self.unit_cell_vertices_count*n
            last=first+self.unit_cell_vertices_count
            polygon = [(self.bent_xs[i], self.bent_ys[i]) for i in range(first,last)]
            
            if self.fillet_radius is not None:
                # create path object, fillet path at control vertices (noncollinear points)
                # then convert to an lwpolyline dxf file entry and add to the file
                p = ezdxf.path.from_vertices(polygon)
                p.close()
                pp = ezdxf.path.fillet(p.control_vertices(), radius=self.fillet_radius)
                out = ezdxf.path.to_lwpolylines((pp,))
                self.msp.add_lwpolyline(next(out), close=True)
            else:
                # no fillets: polygon point list goes straight into the dxf file
                self.msp.add_lwpolyline(polygon, close=True)
            
            bar.next()
        bar.finish()

    def write_dxf(self,filename):
        # Save the DXF file
        self.doc.saveas(filename)

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print('usage: python bender.py /path/to/config.json')
        raise ValueError()

    cfgdict=None
    dxf_filename=None
    with open(sys.argv[1]) as cfgfile:
        cfgdict = json.load(cfgfile)
        dxf_filename = cfgfile.name.split('.json')[0] + '.dxf'
        print(cfgfile.name, cfgdict)

    cfg_tline, cfg_track = cfgdict['tline'], cfgdict['track']

    floquet = tline.FloquetUnitCell()
    for w_loaded_line, n_fishbone_cell in zip(cfg_tline['w_loaded_lines'], cfg_tline['n_fishbone_cells']):
        fishbone = tline.FishboneUnitCell(cfg_tline['fishbone_cell_length'], cfg_tline['w_center_line'], cfg_tline['w_fishbone_line'], w_loaded_line, cfg_tline['gnd_spacing'], cfg_tline['gnd_interdigitate'])
        floquet.append_fishbones(fishbone, n_fishbone_cell)

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


    trackseq = track.TrackSequence()
    if cfg_track['style'] == 'fermat':
        min_spacing = cfg_track['line_spacing_factor']*floquet.cell_min_spacing()
        turns = cfg_track['n_turns']
        final_angle = turns * 2 * np.pi

        tline_compact_length = cfg_track['compact_length']

    # fermat1 = track.FermatSpiralTrack(turns, min_spacing, True)
    # fermat2 = track.FermatSpiralTrack(turns, min_spacing, False)
    # arc2 = fermat2.construct_suitable_ArcTrack()
    # arc1 = fermat1.construct_suitable_ArcTrack()
    # straight2 = arc2.construct_suitable_StraightTrack(tline_compact_length)
    # straight1 = arc1.construct_suitable_StraightTrack(tline_compact_length)

        trackseq.append_track(straight2)
        trackseq.append_track(arc2)
        trackseq.append_track(fermat2)
        trackseq.append_track(fermat1)
        trackseq.append_track(arc1)
        trackseq.append_track(straight1)
    elif cfg_track['style'] == 'singletary':
        # filleted square wave shape

        # math and variable assignment
        compact_width = cfg_track['compact_width']
        compact_height = cfg_track['compact_height']
        r_launch = cfg_track['radius_launch']
        launch_len = cfg_track['launch_len']
        n_half_periods = cfg_track['n_half_periods']
        # working-width
        working_width = compact_width - (2 * r_launch) - (2 * launch_len)
        # half-period 'T'
        T = working_width / n_half_periods
        # semi-circle radius 'R'
        R = T / 2
        # amplitude of straight-line section
        A = (compact_height / 2) - R

        # Building track
        x_start = launch_len + r_launch
        x_end = x_start + n_half_periods * T

        launchstraight1 = track.StraightTrack(np.array((0,0)), np.array((launch_len,0)))
        launchcurve1 = track.ArcTrack(np.array((launch_len, r_launch)), r_launch, 3*np.pi/2, 2*np.pi, False)
        launchstraight2 = track.StraightTrack(np.array((launch_len+r_launch, r_launch)), np.array((launch_len+r_launch, A)))

        delaunchstraight1 = track.StraightTrack(np.array((x_end, -A)), np.array((x_end, -r_launch)))
        delaunchcurve1 = track.ArcTrack(np.array((x_end+r_launch, -r_launch)), r_launch, -np.pi, -np.pi/2, True)
        delaunchstraight2 = track.StraightTrack(np.array((x_end + r_launch, 0)),
                                              np.array((x_end + r_launch + launch_len, 0)))

        trackseq.append_track(launchstraight1)
        trackseq.append_track(launchcurve1)
        trackseq.append_track(launchstraight2)

        if n_half_periods % 2 != 0:
            print('n_half_periods must be even')
            ValueError()

        for i in range(n_half_periods):
            curve, straight = None, None
            if i % 2 == 0:
                curve = track.ArcTrack(np.array((x_start+(0.5 + i) * T, A)), R, -np.pi, 0, True)
                straight = track.StraightTrack(np.array((x_start+(i+1)*T,A)), np.array((x_start+(i+1)*T,-A)))
            else:
                curve = track.ArcTrack(np.array((x_start+(0.5 + i) * T, -A)), R, np.pi, 2*np.pi, False)
                straight = track.StraightTrack(np.array((x_start+(i + 1) * T, -A)), np.array((x_start+(i + 1) * T, A)))
            trackseq.append_track(curve)
            if i < n_half_periods-1:
                trackseq.append_track(straight)

        trackseq.append_track(delaunchstraight1)
        trackseq.append_track(delaunchcurve1)
        trackseq.append_track(delaunchstraight2)

    bender = Bender(trackseq, floquet, cfg_tline['fillet_radius'])
    print('construct_tline')
    bender.construct_tline()
    print('bend_tline')
    bender.bend_tline(plot=False)
    print('write_bent_tline')
    bender.write_bent_tline()

    print('mirror_tline')
    bender.mirror_tline()
    print('bend_tline')
    bender.bend_tline(plot=False)
    print('write_bent_tline')
    bender.write_bent_tline()

    print(f'{floquet.cell_length()*1e6:.1f}um unit cell * {bender.n_floquet_cells:d} unit cells should = {trackseq.total_arclength()*1e6:.0f}um track length')

    print('write_dxf')
    bender.write_dxf(dxf_filename)

