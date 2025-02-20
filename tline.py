# load .dat file with unit cell parameters
# construct polygon arrays for each relevant line
# break into fishbone unit cells

# construct spiral track polygon (Fermat spiral with diametrically opposed launch curves)
#     parametrized by arclength
# sample spiral track at arclengths given by x-values of artwork polygon fishbone cell starts
# compute spiral track tangents and normals at same sample arclengths

# offset from track sample points by fishbone unit cell rotated into (tangent, normal) coordinate system

# generate svg file 

import numpy as np
import matplotlib.pyplot as plt
import ezdxf
import ezdxf.path
from ezdxf.math import Vec3


class FishboneUnitCell:
    def __init__(self, cell_length, w_center_line, w_fishbone_line, w_line):
        self.cell_length, self.w_center_line, self.w_fishbone_line, self.w_line = cell_length, w_center_line, w_fishbone_line, w_line

    def vertices(self):
        xvals1 = [0,
                  self.cell_length / 2 - self.w_line / 2, self.cell_length / 2 - self.w_line / 2,
                  self.cell_length / 2 + self.w_line / 2, self.cell_length / 2 + self.w_line / 2,
                  self.cell_length]#, self.cell_length,
                  # self.cell_length / 2 + self.w_line / 2, self.cell_length / 2 + self.w_line / 2,
                  # self.cell_length / 2 - self.w_line / 2, self.cell_length / 2 - self.w_line / 2,
                  # 0]

        yvals1 = [self.w_center_line / 2,
                  self.w_center_line / 2, self.w_fishbone_line + self.w_center_line / 2,
                  self.w_fishbone_line + self.w_center_line / 2, self.w_center_line / 2,
                  self.w_center_line / 2]#, -self.w_center_line / 2,
                  # -self.w_center_line / 2, -self.w_fishbone_line - self.w_center_line / 2,
                  # -self.w_fishbone_line - self.w_center_line / 2, -self.w_center_line / 2,
                  # -self.w_center_line / 2
                  # ]

        return np.array((xvals1, yvals1))


class FloquetUnitCell:
    def __init__(self):
        self.fishbones = []

    def append_fishbones(self, fishbone_cell, n_fishbones=1):
        for n in range(n_fishbones):
            self.fishbones = np.append(self.fishbones, fishbone_cell)

    def vertices(self):
        v = np.array([[], []])
        xstart = 0
        for fishbone in self.fishbones:
            fishbone_vertices = fishbone.vertices()
            fishbone_vertices[0, :] += xstart
            v = np.append(v, fishbone_vertices, axis=1)
            xstart += fishbone.cell_length

        for fishbone in reversed(self.fishbones):
            xstart -= fishbone.cell_length
            reverse_vertices = fishbone.vertices()
            reverse_vertices = reverse_vertices[:, ::-1]
            reverse_vertices[1] = -reverse_vertices[1]
            reverse_vertices[0,:]+=xstart
            v=np.append(v, reverse_vertices, axis=1)

        v = np.append(v, v[:,0:1], axis=1)
        return v

    def cell_length(self):
        return sum([fishbone.cell_length for fishbone in self.fishbones])

    def cell_min_spacing(self):
        return max(np.concatenate([fishbone.vertices()[1, :] for fishbone in self.fishbones]))


if __name__ == '__main__':
    pass
    import track

    fishboneA = FishboneUnitCell(3e-6, 1.5e-6, 42e-6, 1.5e-6, gnd_spacing=10e-6, interdigitate=False)
    ustrip = FloquetUnitCell()
    ustrip.append_fishbones(fishboneA, 14)
    xs, ys = ustrip.vertices()

    straight1 = track.StraightTrack(np.array([0, 0]), np.array([1e-6, 0]))

    fig, ax = plt.subplots()
    ax.plot(xs, ys)
    # for i in range(len(xs)):
    #     ax.annotate(str(i), (xs[i],ys[i]))
    plt.show()

    merged_list = [(xs[i], ys[i]) for i in range(0, len(xs))]
    # merged_list_gnd = [(xgnd[i], ygnd[i]) for i in range(0, len(xgnd))]

    doc = ezdxf.new()
    msp = doc.modelspace()

    p = ezdxf.path.from_vertices(merged_list)
    p.close()
    pp = ezdxf.path.fillet(p.control_vertices(), radius=0.75e-6)

    out = ezdxf.path.to_lwpolylines({pp})
    for n in out:
        msp.add_lwpolyline(n, close=True)

    # out = ezdxf.path.render_lwpolylines(layout=msp, paths={pp})

    # Save the DXF file
    doc.saveas("unitcell.dxf")

    # print(merged_list)
