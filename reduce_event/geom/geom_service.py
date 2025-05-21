# geom/geom_service.py

import pandas as pd

CHAM_LUT_MAP = {
    31: [1, 2, 3, 4, 5, 6],        # H1B → D0U, D0Up, D0X, D0Xp, D0V, D0Vp
    32: [1, 2, 3, 4, 5, 6],        # H1T → same as above

    37: [13, 14, 15, 16, 17, 18],  # H2B → D2V, D2Vp, D2Xp, D2X, D2U, D2Up
    38: [13, 14, 15, 16, 17, 18],  # H2T → same as above

    40: [19, 20, 21, 22, 23, 24],  # H3T → D3pVp, D3pV, D3pXp, D3pX, D3pUp, D3pU
    39: [25, 26, 27, 28, 29, 30],  # H3B → D3mVp, D3mV, D3mXp, D3mX, D3mUp, D3mU
}

TX_MAX = 0.3  # maximum projected slope in X for hodo-to-chamber mapping

class Plane:
    def __init__(self, detectorID, planeType, z0, n_elements, spacing, cellWidth):
        self.detectorID = detectorID
        self.planeType = planeType
        self.z0 = z0
        self.n_elements = n_elements
        self.spacing = spacing
        self.cellWidth = cellWidth
        self.xoffset = 0.0
        self.x0 = 0.0
        self.y0 = 0.0
        self.costheta = 1.0
        self.sintheta = 0.0

    def get_wire_position(self, elementID):
        return (elementID - (self.n_elements + 1) / 2.0) * self.spacing + self.xoffset

    def get_2d_box_size(self, elementID):
        x_center = self.get_wire_position(elementID)
        return (x_center - 0.5 * self.cellWidth,
                x_center + 0.5 * self.cellWidth,
                -50.0, 50.0)


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
                planeType=1,
                z0=row.z0,
                n_elements=row.n_ele,
                spacing=row.cell_spacing,
                cellWidth=row.cell_width
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
        Constructs a lookup table mapping chamber hits to hodos that can mask them.

        - For each hodoscope paddle, we project its x-extent (box size) to every chamber z-plane.
        - If a chamber wire center falls inside that projected x-range, we record the mapping.

        NOTE:
        - This is chamber-to-hodo (c2h) — used to decide whether to keep a chamber hit based on hodo hits.
        - If a hodo ID or chamber ID is missing, we skip it with a warning.
        """
        for hodo_id, cham_ids in CHAM_LUT_MAP.items():
            if hodo_id not in self.detectors:
                print(f"[WARNING] Hodo {hodo_id} not in geometry — skipping.")
                continue

            z_hodo = self.get_plane_position(hodo_id)

            for hodo_elem in range(1, self.get_n_elements(hodo_id) + 1):
                hodo_uid = hodo_id * 1000 + hodo_elem
                x0_min, x0_max, _, _ = self.get_2d_box_size(hodo_id, hodo_elem)

                for cham_id in cham_ids:
                    if cham_id not in self.detectors:
                        print(f"[WARNING] Chamber {cham_id} not in geometry — skipping.")
                        continue

                    z_cham = self.get_plane_position(cham_id)
                    dz = z_cham - z_hodo
                    x_min = x0_min - abs(TX_MAX * dz)
                    x_max = x0_max + abs(TX_MAX * dz)

                    for cham_elem in range(1, self.get_n_elements(cham_id) + 1):
                        x_c = self.detectors[cham_id].get_wire_position(cham_elem)
                        if x_min <= x_c <= x_max:
                            cham_uid = cham_id * 1000 + cham_elem
                            self.c2h.setdefault(cham_uid, []).append(hodo_uid)
