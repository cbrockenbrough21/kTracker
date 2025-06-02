# geom/geom_service.py

import pandas as pd
import math

CHAM_LUT_MAP = {
    31: [1, 2, 3, 4, 5, 6],        # H1B → D0U, D0Up, D0X, D0Xp, D0V, D0Vp
    32: [1, 2, 3, 4, 5, 6],        # H1T → same as above

    37: [13, 14, 15, 16, 17, 18],  # H2B → D2V, D2Vp, D2Xp, D2X, D2U, D2Up
    38: [13, 14, 15, 16, 17, 18],  # H2T → same as above

    40: [19, 20, 21, 22, 23, 24],  # H3T → D3pVp, D3pV, D3pXp, D3pX, D3pUp, D3pU
    39: [25, 26, 27, 28, 29, 30],  # H3B → D3mVp, D3mV, D3mXp, D3mX, D3mUp, D3mU
}

TX_MAX = 0.15  # maximum projected slope in X for hodo-to-chamber mapping
TY_MAX = 0.10  # maximum projected slope in Y for hodo-to-chamber mapping
BUFFER = 2  # number of elements to buffer on each side of matched range

def get_plane_type(det_name: str, angle_from_vert: float) -> int:
    if any(k in det_name for k in ["X", "T", "B"]):
        return 1
    elif any(k in det_name for k in ["U", "V"]):
        return 2 if angle_from_vert > 0 else 3
    elif any(k in det_name for k in ["Y", "L", "R"]):
        return 4
    else:
        return -1  # Unknown
    
# def get_plane_type(det_name: str) -> int:
#     if det_name.startswith("D"):
#         return 1
#     elif det_name.startswith("H"):
#         return 0
#     elif det_name.startswith("P"):
#         return 2
#     elif det_name.startswith("DP"):
#         return 3
#     else:
#         return -1

def line_crossing(x1, y1, x2, y2, x3, y3, x4, y4):
    """
    Checks if segment (x1, y1)-(x2, y2) intersects with (x3, y3)-(x4, y4).
    """
    tc = (x1 - x2) * (y3 - y1) + (y1 - y2) * (x1 - x3)
    td = (x1 - x2) * (y4 - y1) + (y1 - y2) * (x1 - x4)
    return tc * td < 0

class Plane:
    def __init__(self, detectorID, planeType, x0, y0, z0, n_elements, spacing, cellWidth,
                 angle_from_vert, xoffset, height, theta_x, theta_y, theta_z, deltaW):

        self.detectorID = detectorID
        self.planeType = planeType
        self.n_elements = n_elements
        self.spacing = spacing
        self.cellWidth = cellWidth
        self.angle_from_vert = angle_from_vert  # radians
        self.xoffset = xoffset
        self.width = spacing * n_elements  
        self.height = height
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.theta_x = theta_x
        self.theta_y = theta_y
        self.theta_z = theta_z

        # Derived geometry variables
        self.costheta = math.cos(angle_from_vert + theta_z)
        self.sintheta = math.sin(angle_from_vert + theta_z)
        self.tantheta = math.tan(angle_from_vert + theta_z)

        self.y1 = y0 - 0.5 * height
        self.y2 = y0 + 0.5 * height
        self.deltaW = deltaW  
        self.deltaX = self.deltaW * self.costheta  
    
    def get_wire_endpoints(self, elementID):
        """
        Returns the (x_min, x_max, y_min, y_max) endpoints of a slanted wire in 2D
        following the same logic as GeomSvc::getWireEndPoints in C++
        """
        y_min = self.y1
        y_max = self.y2

        mid_index = (self.n_elements + 1) / 2.0
        dw = self.xoffset + self.spacing * (elementID - mid_index)

        x_center = (self.x0 + self.deltaX) + dw * math.cos(self.theta_z)
        tan_theta = math.tan(self.theta_z)
        dx = 0.5 * abs(tan_theta * (y_max - y_min))

        x_min = x_center - dx
        x_max = x_center + dx

        return x_min, x_max, y_min, y_max

    def get_wire_position(self, elementID):
        mid_index = (self.n_elements + 1) / 2.0
        dw = self.spacing * (elementID - mid_index) + self.xoffset
        angle = self.theta_z  # already in radians

        x_pos = self.x0 * math.cos(angle) + self.y0 * math.sin(angle) + dw + self.deltaW
        return x_pos

    def get_2d_box_size(self, elementID):
        """
        Python translation of GeomSvc::get2DBoxSize.
        Returns (x_min, x_max, y_min, y_max) for the element's 2D box.
        """
        if self.planeType == 1:
            # Get x_center from wire position
            x_center = self.get_wire_position(elementID)
            x_width = 0.5 * self.cellWidth

            x_min = x_center - x_width
            x_max = x_center + x_width

            y_min = self.y1
            y_max = self.y2
        else:
            # For slanted planes, box is vertical
            y_center = self.get_wire_position(elementID)
            y_width = 0.5 * self.cellWidth

            y_min = y_center - y_width
            y_max = y_center + y_width

            x_min = self.x0 - 0.5 * self.width
            x_max = self.x0 + 0.5 * self.width

        return x_min, x_max, y_min, y_max


class GeometryService:
    def __init__(self, tsv_path: str):
        self.detectors = {}  # detectorID -> Plane instance
        self.c2h = {}        # chamberUID -> list of masking hodoUIDs
        self.tsv_path = tsv_path
        self.load_geometry_from_tsv()

    def load_geometry_from_tsv(self):
        columns = [
            "det_name", "n_ele", "cell_spacing", "cell_width", "angle_from_vert",
            "xoffset", "width", "height", "x0", "y0", "z0", "theta_x", "theta_y", "theta_z"
        ]
        df = pd.read_csv(self.tsv_path, sep='\t', comment='#', names=columns)

        for det_id, row in enumerate(df.itertuples(), start=1):
            plane = Plane(
                detectorID=det_id,
                planeType=get_plane_type(row.det_name, row.angle_from_vert),
                x0=row.x0,
                y0=row.y0,
                z0=row.z0,
                height=row.height,
                n_elements=row.n_ele,
                spacing=row.cell_spacing,
                cellWidth=row.cell_width,
                angle_from_vert=row.angle_from_vert,
                xoffset=row.xoffset,
                theta_x=row.theta_x,
                theta_y=row.theta_y,
                theta_z=row.theta_z,
                deltaW=0.0   # Assuming deltaW is not provided in the TSV, set to 0.0
            )
            self.detectors[det_id] = plane

        self.init_hodo_mask_lut()

    def get_plane_position(self, det_id):
        return self.detectors[det_id].z0

    def get_plane_type(self, det_id):
        return self.detectors[det_id].planeType

    def get_n_elements(self, det_id):
        return self.detectors[det_id].n_elements

    def get_2d_box_size(self, det_id, elem_id):
        return self.detectors[det_id].get_2d_box_size(elem_id)

    def init_hodo_mask_lut(self):
        """
        Constructs a lookup table mapping chamber hits (by UID) to hodoscope hits that can justify keeping them.
        Includes TX_MAX and TY_MAX projections and ±BUFFER element padding.
        """
        for hodo_id, cham_ids in CHAM_LUT_MAP.items():
            if hodo_id not in self.detectors:
                print(f"[WARNING] Hodo {hodo_id} not in geometry — skipping.")
                continue

            z_hodo = self.get_plane_position(hodo_id)

            for hodo_elem in range(1, self.get_n_elements(hodo_id) + 1):
                hodo_uid = hodo_id * 1000 + hodo_elem
                x0_min, x0_max, y0_min, y0_max = self.get_2d_box_size(hodo_id, hodo_elem)

                for cham_id in cham_ids:
                    if cham_id not in self.detectors:
                        print(f"[WARNING] Chamber {cham_id} not in geometry — skipping.")
                        continue

                    z_cham = self.get_plane_position(cham_id)
                    dz = z_cham - z_hodo

                    # Expand hodo paddle projection to chamber Z using TX_MAX and TY_MAX margin
                    x_min = x0_min - abs(TX_MAX * dz)
                    x_max = x0_max + abs(TX_MAX * dz)
                    y_min = y0_min - abs(TY_MAX * dz)
                    y_max = y0_max + abs(TY_MAX * dz)

                    n_elements = self.get_n_elements(cham_id)

                    # Store (eid, x_pos, y_min, y_max) for chamber wires
                    wire_info = [
                        (
                            eid,
                            self.detectors[cham_id].get_wire_position(eid),
                            self.detectors[cham_id].y0 - 0.5 * self.detectors[cham_id].height,
                            self.detectors[cham_id].y0 + 0.5 * self.detectors[cham_id].height,
                        )
                        for eid in range(1, n_elements + 1)
                    ]
                
                    # Find element range where x and y projections overlap
                    lo = hi = None
                    for idx, (eid, x, y_lo, y_hi) in enumerate(wire_info):
                        if x_min <= x <= x_max and y_min <= y_hi and y_max >= y_lo:
                            if lo is None:
                                lo = idx
                            hi = idx

                    if lo is None or hi is None:
                        continue

                    # Apply ±2 element buffer
                    lo = max(0, lo - BUFFER)
                    hi = min(n_elements - 1, hi + BUFFER)

                    # Map each chamber UID in this buffered range to the hodo UID
                    for idx in range(lo, hi + 1):
                        cham_elem = wire_info[idx][0]
                        cham_uid = cham_id * 1000 + cham_elem
                        self.c2h.setdefault(cham_uid, []).append(hodo_uid)