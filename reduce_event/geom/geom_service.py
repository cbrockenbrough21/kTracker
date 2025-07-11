# geom/geom_service.py

import pandas as pd
import math
import bisect
from collections import defaultdict
from reduce_event.reco_constants import N_CHAMBER_PLANES, N_HODO_PLANES, N_PROP_PLANES, TX_MAX, TY_MAX, BUFFER

CHAM_LUT_MAP = {
    31: [1, 2, 3, 4, 5, 6],        # H1B → D0U, D0Up, D0X, D0Xp, D0V, D0Vp
    32: [1, 2, 3, 4, 5, 6],        # H1T → same as above

    37: [13, 14, 15, 16, 17, 18],  # H2B → D2V, D2Vp, D2Xp, D2X, D2U, D2Up
    38: [13, 14, 15, 16, 17, 18],  # H2T → same as above

    40: [19, 20, 21, 22, 23, 24],  # H3T → D3pVp, D3pV, D3pXp, D3pX, D3pUp, D3pU
    39: [25, 26, 27, 28, 29, 30],  # H3B → D3mVp, D3mV, D3mXp, D3mX, D3mUp, D3mU
    
    46: [],  # H4B
    47: [],  # H4T
}

def to_local_detector_name(detector_name: str, element_id: int):
    if detector_name.startswith("P") and len(detector_name) > 3 and detector_name[2] in {'H', 'V'}:
        view = detector_name[2]
        fb = detector_name[3:]

        # Determine XY and 1/2
        xy = 'Y' if view == 'H' else 'X'
        side = '1' if 'f' in fb else '2'

        # Extract moduleID (the number right after P1 or P2)
        digits = ''.join(c for c in fb if c.isdigit())
        module_id = int(digits)

        # Build new detector name
        base = detector_name[:2]  # e.g., 'P1'
        new_name = f"{base}{xy}{side}"

        # Adjust element ID if ≤ 8
        if element_id <= 8:
            element_id = (9 - module_id) * 8 + element_id

        return new_name, element_id
    return detector_name, element_id

def line_crossing(x1, y1, x2, y2, x3, y3, x4, y4):
    """
    Python translation of EventReducer::lineCrossing in C++ (EventReducer.cxx).
    Checks if segment (x1, y1)-(x2, y2) intersects with (x3, y3)-(x4, y4).
    """
    tc = (x1 - x2) * (y3 - y1) + (y1 - y2) * (x1 - x3)
    td = (x1 - x2) * (y4 - y1) + (y1 - y2) * (x1 - x4)
    return tc * td < 0

class Plane:
    def __init__(self, detector_id, detector_name, x0, y0, z0, n_elements, spacing, cell_width,
                 angle_from_vert, xoffset, height, theta_x, theta_y, theta_z, delta_w,
                 plane_type=None):  # Plane type is optional
        # --- Detector identifier ---
        self.detector_id = detector_id
        self.detector_name = detector_name

        # --- Ideal properties ---
        self.n_elements = n_elements
        self.spacing = spacing
        self.cell_width = cell_width
        self.xoffset = xoffset
        self.overlap = cell_width - spacing
        self.angle_from_vert = angle_from_vert  # radians
        
        # Infer planeType if not provided
        if plane_type is None:
            self.planeType = self.infer_plane_type()
        else:
            self.planeType = plane_type

        # Trigonometric setup using rZ (future-proofed for alignment)
        self.rZ = theta_z  # rotZ defaults to 0 for now
        self.costheta = math.cos(angle_from_vert + self.rZ)
        self.sintheta = math.sin(angle_from_vert + self.rZ)
        self.tantheta = math.tan(angle_from_vert + self.rZ)

        # --- Survey info ---
        self.x0 = x0  # center of detector
        self.y0 = y0
        self.z0 = z0

        self.width = spacing * n_elements
        self.height = height

        self.x1 = self.x0 - 0.5 * self.width
        self.x2 = self.x0 + 0.5 * self.width
        self.y1 = self.y0 - 0.5 * height
        self.y2 = self.y0 + 0.5 * height
        self.z1 = None  # optional
        self.z2 = None  # optional

        self.theta_x = theta_x
        self.theta_y = theta_y
        self.theta_z = theta_z

        # --- Alignment info ---
        self.delta_w = delta_w
        self.deltaX = self.delta_w * self.costheta
        self.deltaY = self.delta_w * self.sintheta
        self.deltaZ = 0.0
        self.deltaW_module = [0.0] * 9  # for prop tubes, unused
        self.rotX = 0.0
        self.rotY = 0.0
        self.rotZ = 0.0
        self.resolution = 0.0

        # --- Final position/rotation (unused now) ---
        self.xc = self.x0 + self.deltaX
        self.yc = self.y0 + self.deltaY
        self.zc = self.z0 + self.deltaZ
        self.wc = self.xc * self.costheta + self.yc * self.sintheta
        self.rX = self.theta_x + self.rotX
        self.rY = self.theta_y + self.rotY
        self.rZ = self.theta_z + self.rotZ  # already set above, redundancy ok

        # --- Geometric setup (unused vectors for now) ---
        self.nVec = None
        self.uVec = None
        self.vVec = None
        self.xVec = None
        self.yVec = None
        self.rotM = None

        # --- Calibration info (not used yet) ---
        self.tmin = None
        self.tmax = None
        self.rtprofile = None

        # --- Wire positions (like C++ elementPos) ---
        self.elementPos = [
            ((i - (n_elements + 1) / 2.0) * spacing + xoffset +
             x0 * self.costheta + y0 * self.sintheta + delta_w)
            for i in range(1, n_elements + 1)
        ]
        #self.elementPos.sort() - doesnt change accuracy 
    
    def infer_plane_type(self):
        """
        GeomSvc.cxx lines 481-497
        """
        if any(k in self.detector_name for k in ["X", "T", "B"]):
            return 1
        elif any(k in self.detector_name for k in ["U", "V"]):
            return 2 if self.angle_from_vert > 0 else 3
        elif any(k in self.detector_name for k in ["Y", "L", "R"]):
            return 4
        else:
            return -1
    
    def get_wire_endpoints(self, elementID):
        """
        Python translation of GeomSvc::getWireEndPoints() in C++ (GeomSvc.cxx).
        
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
        """
        Python translation of GeomSvc::getWirePosition() in C++ (GeomSvc.cxx).

        Calculates the projected wire position in the rotated frame using
        the detector center (x0, y0), slant angle, and alignment offset deltaW.
        """
        
        mid_index = (self.n_elements + 1) / 2.0
        dw = (elementID - mid_index) * self.spacing + self.xoffset #displacement of wire
        x_pos = dw + self.x0 * self.costheta + self.y0 * self.sintheta + self.delta_w
        return x_pos

    def get_2d_box_size(self, elementID):
        """
        Python translation of GeomSvc::get2DBoxSize.
        Returns (x_min, x_max, y_min, y_max) for the element's 2D box.
        """
        if self.planeType == 1:
            # Get x_center from wire position
            x_center = self.get_wire_position(elementID)
            x_width = 0.5 * self.cell_width

            x_min = x_center - x_width
            x_max = x_center + x_width

            y_min = self.y1
            y_max = self.y2
        else:
            # For slanted planes, box is vertical
            y_center = self.get_wire_position(elementID)
            y_width = 0.5 * self.cell_width

            y_min = y_center - y_width
            y_max = y_center + y_width

            x_min = self.x0 - 0.5 * self.width
            x_max = self.x0 + 0.5 * self.width

        return x_min, x_max, y_min, y_max


class GeometryService:
    def __init__(self, tsv_path: str):
        self.detectors = {}  # detectorID -> Plane instance
        self.c2h = {}        # chamberUID -> list of masking hodoUIDs
        self.h2celementID_lo = defaultdict(list)
        self.h2celementID_hi = defaultdict(list)
        self.c2helementIDs = defaultdict(list)
        self.CHAM_LUT_MAP = {
            31: [1, 2, 3, 4, 5, 6],
            32: [1, 2, 3, 4, 5, 6],
            37: [13, 14, 15, 16, 17, 18],
            38: [13, 14, 15, 16, 17, 18],
            40: [19, 20, 21, 22, 23, 24],
            39: [25, 26, 27, 28, 29, 30],
            46: [],
            47: [],
        }
        self.tsv_path = tsv_path
        self.load_geometry_from_tsv()

    def load_geometry_from_tsv(self):
        columns = [
            "det_name", "n_ele", "cell_spacing", "cell_width", "angle_from_vert",
            "xoffset", "width", "height", "x0", "y0", "z0", "theta_x", "theta_y", "theta_z"
        ]
        df = pd.read_csv(self.tsv_path, sep='\t', comment='#', names=columns)
        
        name_to_id = {
            "D0U": 1, "D0Up": 2, "D0X": 3, "D0Xp": 4, "D0V": 5, "D0Vp": 6,
            "D1V": 7, "D1Vp": 8, "D1X": 9, "D1Xp": 10, "D1U": 11, "D1Up": 12,
            "D2V": 13, "D2Vp": 14, "D2Xp": 15, "D2X": 16, "D2U": 17, "D2Up": 18,
            "D3pVp": 19, "D3pV": 20, "D3pXp": 21, "D3pX": 22, "D3pUp": 23, "D3pU": 24,
            "D3mVp": 25, "D3mV": 26, "D3mXp": 27, "D3mX": 28, "D3mUp": 29, "D3mU": 30,
            "H1B": 31, "H1T": 32, "H1L": 33, "H1R": 34, "H2L": 35, "H2R": 36,
            "H2B": 37, "H2T": 38, "H3B": 39, "H3T": 40,
            "H4Y1L": 41, "H4Y1R": 42, "H4Y2L": 43, "H4Y2R": 44, "H4B": 45, "H4T": 46,
            "P1Y1": 47, "P1Y2": 48, "P1X1": 49, "P1X2": 50,
            "P2X1": 51, "P2X2": 52, "P2Y1": 53, "P2Y2": 54,
            "DP1TL": 55, "DP1TR": 56, "DP1BL": 57, "DP1BR": 58,
            "DP2TL": 59, "DP2TR": 60
        }

        for row in df.itertuples():
            det_name = row.det_name.strip()
            if det_name not in name_to_id:
                #print(f"[WARNING] Detector name '{det_name}' not in name-to-ID mapping. Skipping.")
                continue
            else:
                det_id = name_to_id[det_name]
            plane = Plane(
                detector_id=det_id,
                detector_name=row.det_name,
                x0=row.x0,
                y0=row.y0,
                z0=row.z0,
                height=row.height,
                n_elements=row.n_ele,
                spacing=row.cell_spacing,
                cell_width=row.cell_width,
                angle_from_vert=row.angle_from_vert,
                xoffset=row.xoffset,
                theta_x=row.theta_x,
                theta_y=row.theta_y,
                theta_z=row.theta_z,
                delta_w=0.0   # Assuming deltaW is not provided in the TSV, set to 0.0 
            )
            self.detectors[det_id] = plane

        self.init_hodo_mask_lut()
        #self.init_hodo_mask_lut_new()
    
    def dump_geometry_summary(self, output_path="geometry_dump.tsv"):
        rows = []
        for det_id, plane in self.detectors.items():
            row = {
                "detector_id": det_id,
                "detector_name": plane.detector_name,
                "plane": plane.planeType,
                "n_elements": plane.n_elements,
                "spacing": plane.spacing,
                "xoffset": plane.xoffset,
                "x0": plane.x0,
                "y0": plane.y0,
                "z0": plane.z0,
                "angle_from_vert": plane.angle_from_vert,
                "elementPos_1": plane.elementPos[0] if plane.elementPos else None,
                "elementPos_2": plane.elementPos[1] if len(plane.elementPos) > 1 else None,
                "elementPos_last": plane.elementPos[-1] if plane.elementPos else None,
            }
            rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(output_path, sep="\t", index=False)
        print(f"[INFO] Geometry summary dumped to {output_path}")       

    def get_plane_position(self, det_id):
        return self.detectors[det_id].z0

    def get_plane_type(self, det_id):
        return self.detectors[det_id].planeType

    def get_plane_n_elements(self, det_id):
        return self.detectors[det_id].n_elements

    def get_2d_box_size(self, det_id, elem_id):
        return self.detectors[det_id].get_2d_box_size(elem_id)
    
    def get_wire_endpoints(self, det_id, elem_id):
        return self.detectors[det_id].get_wire_endpoints(elem_id)
    
    
    
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

            for hodo_elem in range(1, self.get_plane_n_elements(hodo_id) + 1):
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

                    n_elements = self.get_plane_n_elements(cham_id)

                    plane_type = self.get_plane_type(cham_id)

                    if plane_type == 1:
                        elementID_lo = self.get_exp_element_id(cham_id, x_min)
                        elementID_hi = self.get_exp_element_id(cham_id, x_max)
                                               
                    else:
                        if cham_id in [13, 14, 15, 16, 17, 18]:
                            print(f"\n== Debug cham_id: {cham_id} ==")
                            print(f"x_min = {x_min}, x_max = {x_max}")
                            print(f"y_min = {y_min}, y_max = {y_max}")

                        detector = self.detectors[cham_id]

                        # Run both versions
                        el_lo_len, el_hi_len = get_element_range_lenient(detector, n_elements, x_min, x_max, y_min, y_max)
                        el_lo_cross, el_hi_cross = get_element_range_crossing(detector, cham_id, n_elements, x_min, x_max, y_min, y_max)

                        # Print side-by-side debug output
                        if cham_id in [13, 14, 15, 16, 17, 18]:
                            print(f"[Lenient] → Range: {el_lo_len}–{el_hi_len}" if el_lo_len else "[Lenient] → No elements found")
                            print(f"[Crossing] → Range: {el_lo_cross}–{el_hi_cross}" if el_lo_cross else "[Crossing] → No elements found")
                        
                        # elementID_lo = el_lo_len
                        # elementID_hi = el_hi_len
                        
                        # if elementID_lo is None or elementID_hi is None:
                        #     continue
                        
                        if el_lo_len is not None and el_lo_cross is not None:
                            elementID_lo = min(el_lo_len, el_lo_cross)
                        elif el_lo_len is not None:
                            elementID_lo = el_lo_len
                        elif el_lo_cross is not None:
                            elementID_lo = el_lo_cross
                        else:
                            continue  # Skip this chamber if no lower bound

                        if el_hi_len is not None and el_hi_cross is not None:
                            elementID_hi = max(el_hi_len, el_hi_cross)
                        elif el_hi_len is not None:
                            elementID_hi = el_hi_len
                        elif el_hi_cross is not None:
                            elementID_hi = el_hi_cross
                        else:
                            continue  # Skip this chamber if no upper bound

                    # Apply ±BUFFER
                    elementID_lo = max(1, elementID_lo - BUFFER)
                    elementID_hi = min(n_elements, elementID_hi + BUFFER)

                    # Map range to hodo UID
                    for eid in range(elementID_lo, elementID_hi + 1):
                        cham_uid = cham_id * 1000 + eid
                        self.c2h.setdefault(cham_uid, []).append(hodo_uid)

    def init_hodo_mask_lut_new(self):
        """
        Python translation of EventReducer::initHodoMaskLUT()

        Initializes hodoscope-to-chamber and chamber-to-hodoscope LUTs:
            - self.h2celementID_lo
            - self.h2celementID_hi
            - self.c2helementIDs
        """
        self.h2celementID_lo.clear()
        self.h2celementID_hi.clear()
        
        for hodo_id, cham_ids in self.CHAM_LUT_MAP.items():
            if hodo_id not in self.detectors:
                print(f"[SKIP] Hodo {hodo_id} not in geometry.")
                continue

            print(f"[INFO] Processing hodo {hodo_id} with {len(cham_ids)} chambers.")
            
            n_paddles = self.get_plane_n_elements(hodo_id)
            for paddle_id in range(1, n_paddles + 1):
                hodo_uid = hodo_id * 1000 + paddle_id

                z0 = self.get_plane_position(hodo_id)
                x0_min, x0_max, y0_min, y0_max = self.get_2d_box_size(hodo_id, paddle_id)

                for cham_id in cham_ids:
                    if cham_id not in self.detectors:
                        print(f"[WARN] Chamber {cham_id} not in geometry.")
                        continue
                    if cham_id < 1:
                        continue
                    
                    plane = self.detectors[cham_id]
                    #print(f"[DEBUG] Detector {cham_id}: elementPos[:3] = {plane.elementPos[:3]}, ... elementPos[-3:] = {plane.elementPos[-3:]}")

                    z = self.get_plane_position(cham_id)
                    x_min = x0_min - abs(TX_MAX * (z - z0))
                    x_max = x0_max + abs(TX_MAX * (z - z0))
                    y_min = y0_min - abs(TY_MAX * (z - z0))
                    y_max = y0_max + abs(TY_MAX * (z - z0))

                    elementID_lo = self.get_plane_n_elements(cham_id)
                    elementID_hi = 0
                    
                    n_elements = self.get_plane_n_elements(cham_id)

                    if self.get_plane_type(cham_id) == 1:
                        elementID_lo = self.get_exp_element_id(cham_id, x_min)
                        elementID_hi = self.get_exp_element_id(cham_id, x_max)
                        #print(f"[INFO] Straight chamber {cham_id}: x_min={x_min:.2f}, x_max={x_max:.2f} → el_lo={elementID_lo}, el_hi={elementID_hi}")

                    else:
                        if cham_id in [13, 14, 15, 16, 17, 18]:
                            print(f"\n== Debug cham_id: {cham_id} ==")
                            print(f"x_min = {x_min}, x_max = {x_max}")
                            print(f"y_min = {y_min}, y_max = {y_max}")

                        detector = self.detectors[cham_id]

                        # Run both versions
                        el_lo_len, el_hi_len = get_element_range_lenient(detector, n_elements, x_min, x_max, y_min, y_max)
                        el_lo_cross, el_hi_cross = get_element_range_crossing(detector, cham_id, n_elements, x_min, x_max, y_min, y_max)

                        # Print side-by-side debug output
                        if cham_id in [13, 14, 15, 16, 17, 18]:
                            print(f"[Lenient] → Range: {el_lo_len}–{el_hi_len}" if el_lo_len else "[Lenient] → No elements found")
                            print(f"[Crossing] → Range: {el_lo_cross}–{el_hi_cross}" if el_lo_cross else "[Crossing] → No elements found")
                        
                        elementID_lo = el_lo_len
                        elementID_hi = el_hi_len
       
                    # Apply ±BUFFER
                    elementID_lo = max(1, elementID_lo - 2)
                    elementID_hi = min(self.get_plane_n_elements(cham_id), elementID_hi + 2)

                    self.h2celementID_lo[hodo_uid].append(cham_id * 1000 + elementID_lo)
                    self.h2celementID_hi[hodo_uid].append(cham_id * 1000 + elementID_hi)

        self.c2helementIDs.clear()
        for hodo_uid, lo_list in self.h2celementID_lo.items():
            for i, cham_uid_lo in enumerate(lo_list):
                cham_uid_hi = self.h2celementID_hi[hodo_uid][i]
                for j in range(cham_uid_lo, cham_uid_hi + 1):
                    self.c2helementIDs[j].append(hodo_uid)
                    
    def get_exp_element_id(self, detector_id, pos_exp):
        """
        Python translation of GeomSvc::getExpElementID in C++ (GeomSvc.cxx).
        """
        plane = self.detectors[detector_id]
        element_pos = plane.elementPos
        if not element_pos:
            return -1  # or raise exception

        pos_min = element_pos[0] - 0.5 * plane.cell_width
        pos_max = element_pos[-1] + 0.5 * plane.cell_width
        
        #print(f"[DEBUG] get_exp_element_id for det {detector_id}: pos_exp={pos_exp:.2f}, pos_min={pos_min:.2f}, pos_max={pos_max:.2f}")

        if pos_exp > pos_max:
            return plane.n_elements + 1
        if pos_exp < pos_min:
            return 0

        index = bisect.bisect_left(element_pos, pos_exp)
        if index < len(element_pos):
            adjustment = -1 if (element_pos[index] - pos_exp) > 0.5 * plane.spacing else 0
            element_id = index + 1 + adjustment
        else:
            element_id = plane.n_elements

        if detector_id > N_CHAMBER_PLANES + N_HODO_PLANES + N_PROP_PLANES:
            bottom = ((detector_id - 55) & 2) > 0
            if bottom:
                element_id = plane.n_elements + 1 - element_id

        #print(f"[DEBUG] get_exp_element_id for det {detector_id}: pos_exp={pos_exp:.2f}, pos_min={pos_min:.2f}, pos_max={pos_max:.2f}, element_id={element_id}")
        #print(f"[DEBUG] elementPos for det {detector_id}: {plane.elementPos[:3]} ... {plane.elementPos[-3:]}")

        return element_id

def get_element_range_lenient(detector, n_elements, x_min, x_max, y_min, y_max):
    wire_info = [
        (
            eid,
            detector.get_wire_position(eid),
            detector.y0 - 0.5 * detector.height,
            detector.y0 + 0.5 * detector.height,
        )
        for eid in range(1, n_elements + 1)
    ]

    elementID_lo = elementID_hi = None
    for eid, x, y_lo, y_hi in wire_info:
        if x_min <= x <= x_max and y_min <= y_hi and y_max >= y_lo:
            if elementID_lo is None:
                elementID_lo = eid
            elementID_hi = eid
    return elementID_lo, elementID_hi

def get_element_range_crossing(detector, cham_id, n_elements, x_min, x_max, y_min, y_max):
    elementID_lo = n_elements + 1
    elementID_hi = 0
    for m in range(1, n_elements + 1):
        x1, x2, y1, y2 = detector.get_wire_endpoints(m)

        crosses_left = line_crossing(x_min, y_min, x_min, y_max, x1, y1, x2, y2)
        crosses_right = line_crossing(x_max, y_min, x_max, y_max, x1, y1, x2, y2)
        
        if cham_id == 18 and (m >= 41 and m<=55):
            print(f"[DEBUG] Element {m}")
            print(f"  Wire: ({x1:.3f}, {y1:.3f}) → ({x2:.3f}, {y2:.3f})")
            print(f"  Crosses left? {crosses_left}")
            print(f"  Crosses right? {crosses_right}")
            print()

        if not crosses_left and not crosses_right:
            continue

        if m < elementID_lo:
            elementID_lo = m
        if m > elementID_hi:
            elementID_hi = m

    if elementID_lo > elementID_hi:
        return None, None
    return elementID_lo, elementID_hi

