"""Lookup table with inverse-distance weighted interpolation for initial guesses

Field-agnostic: works for any combination of key dimensions and value fields.
Loads read-only seed data shipped with the package.
Auto-scales key dimensions by range for distance calculation.
"""

# --------------------------------------------------
# load necessary modules
# --------------------------------------------------
from __future__ import annotations

import json
from importlib.resources import files

import numpy as np
from numpy.typing import NDArray

# --------------------------------------------------
# package data path — tables shipped with the wheel (read-only)
# --------------------------------------------------
_PACKAGE_DATA = files("simbl") / "data"


# --------------------------------------------------
# LookupTable class: initial guess determination
# --------------------------------------------------
class LookupTable:
    """Lookup table with inverse-distance weighted interpolation

    Each instance handles one of the following (model, wall_bc) combinations:
    - Falkner-Skan with adiabatic wall
    - Falkner-Skan with isothermal wall
    - Falkner-Skan-Cooke with adiabatic wall
    - Falkner-Skan-Cooke with isothermal wall

    Seed data is shipped with the package and is read-only.

    Parameters
    ----------
    key_fields : list[str]
        Field names that form the lookup key (e.g. ["mach", "beta"])
    value_fields : list[str]
        Field names for stored values (e.g. ["fpp_wall", "g_wall"])
    default_values : list[float]
        Fallback values when table is empty (e.g. [0.615, 3.69] for adiabatic flat plate)
    fname : str
        JSON filename for lookup table data (e.g. "fsc_lookup.json")
    """

    # --------------------------------------------------
    # initialization of LookupTable instance
    # --------------------------------------------------
    def __init__(
        self,
        key_fields: list[str],
        value_fields: list[str],
        default_values: list[float],
        fname: str,
    ):
        self.key_fields = key_fields
        self.value_fields = value_fields
        self.default_values = np.array(default_values, dtype=np.float64)
        self.fname = fname

        # internal storage
        # _entries: raw dicts from JSON
        # _points: shape (n_entries, n_key_fields) for distance computation
        # _values: shape (n_entries, n_value_fields) for interpolation
        # _scale:  shape (n_key_fields,) — 1/range per dimension, unit if range is zero
        self._entries: list[dict] = []
        self._points: NDArray[np.float64] = np.empty((0, len(key_fields)))
        self._values: NDArray[np.float64] = np.empty((0, len(value_fields)))
        self._scale: NDArray[np.float64] = np.ones(len(key_fields))

        self._load()

    # --------------------------------------------------
    # magic method for len()
    #
    # returns number of entries in the table when using the len() command on a LookupTable instance
    # --------------------------------------------------
    def __len__(self) -> int:
        return len(self._entries)

    # --------------------------------------------------
    # magic method for repr()
    #
    # returns file name and number of entries when using the print() command on a LookupTable instance
    # --------------------------------------------------
    def __repr__(self) -> str:
        return f"LookupTable({self.fname!r}, {len(self)} entries)"

    # --------------------------------------------------
    # load table from package data
    # --------------------------------------------------
    def _load(self) -> None:
        """Load lookup table from package data"""
        try:
            # use read_text + json.loads because package data may be inside a zip/wheel
            text = (_PACKAGE_DATA / self.fname).read_text(encoding="utf-8")
            data = json.loads(text)
            self._entries = data.get("entries", [])
        except (FileNotFoundError, TypeError):
            # no package data found, predict() will use default_values instead
            pass

        # build numpy arrays from entries for interpolation
        self._build_arrays()

    # --------------------------------------------------
    # build numpy arrays for interpolation from entries list
    # --------------------------------------------------
    def _build_arrays(self) -> None:
        """Build numpy arrays for interpolation from entries list"""

        # skip if no entries (arrays will be empty, scale will be ones)
        if not self._entries:
            return

        # initialize empty lists to collect points and values
        points = []
        values = []

        # extract points and values from entries based on key_fields and value_fields
        # example: self._entries = [{"mach": 4.0, "beta": 0.0, "fpp_wall": 0.615, "g_wall": 3.69}, ...]
        # key_fields = ["mach", "beta"]
        # value_fields = ["fpp_wall", "g_wall"]
        # points = [[4.0, 0.0], ...]
        # values = [[0.615, 3.69], ...]
        #
        # entries missing any required field (e.g. from an old cache version) are skipped
        required_fields = set(self.key_fields) | set(self.value_fields)
        for e in self._entries:
            if not required_fields.issubset(e.keys()):
                continue
            points.append([e[k] for k in self.key_fields])
            values.append([e[v] for v in self.value_fields])

        # all entries were stale — leave arrays empty, predict() falls back to default_values
        if not points:
            return

        # convert to numpy arrays (double precision for interpolation)
        self._points = np.array(points, dtype=np.float64)
        self._values = np.array(values, dtype=np.float64)

        # auto-scale: normalize each key dimension by its range
        # dimensions with zero range (single value) get unit scaling
        ranges = self._points.max(axis=0) - self._points.min(axis=0)
        safe_ranges = np.where(ranges > 0, ranges, 1.0)
        self._scale = 1.0 / safe_ranges

    # --------------------------------------------------
    # predict initial guess values
    # --------------------------------------------------
    def predict(self, key_values: dict[str, float]) -> NDArray[np.float64]:
        """Predict initial guess values from tabulated data using inverse-distance weighted interpolation

        Parameters
        ----------
        key_values : dict[str, float]
            Key field values as {field_name: value} for each key_field

        Returns
        -------
        NDArray[np.float64]
            Predicted shooting variable values
        """

        # if no usable entries are loaded, return default values
        # note: default_values are stored in the _TABLE_CONFIGS dictionaries for falkner_skan and falkner_skan_cooke (see initial_guess.py)
        if len(self._points) == 0:
            return self.default_values.copy()

        # build query point (unknown point, where to interpolate to) from key fields
        # example: key_values = {"mach": 4.0, "beta": 0.0} --> query_point = [4.0, 0.0]
        # same shape as one row of _points, so the distances can be computed as the difference between _points and query_point
        query_point = np.array([key_values[k] for k in self.key_fields], dtype=np.float64)

        # compute scaled distance from query_point to each table entry
        #
        # using inverse distance weighting (IDW) method by Shepard (1968):
        #   Shepard, D. (1968), "A two-dimensional interpolation function for
        #   irregularly-spaced data", Proc. 23rd ACM National Conference, 517-524.
        #
        # steps:
        # 1. compute difference between query point and each entry's key values, scale by range to normalize dimensions
        # 2. compute Euclidean distance for each entry
        # 3. if any entry is very close (distance < 0.01), return that entry's values directly (exact match)
        # 4. otherwise, pick k nearest neighbors (k <= 4) and compute weights as inverse of distance, normalize weights, and return weighted average of neighbor values
        #

        # get number of entries
        n_entries = len(self._entries)
        # initialize distances array
        distances = np.empty(n_entries)
        for i in range(n_entries):
            # distance from query point to entry i, e.g. [mach_i - mach_q, beta_i - beta_q]
            diff = self._points[i] - query_point
            # normalize by range (scale) to give equal weight to each dimension regardless of units or range
            scaled_diff = diff * self._scale
            # compute Euclidean distance as sqrt((mach_i - mach_q)^2 * scale_mach^2 + (beta_i - beta_q)^2 * scale_beta^2)
            distances[i] = np.sqrt(np.sum(scaled_diff**2))

        # get minimum value across all scaled distances (minimum distance from query_point to any entry in the table)
        min_dist = np.min(distances)

        # exact match: if any entry is very close to the query point (skip interpolation)
        if min_dist < 0.01:
            # get index of smallest distance (closest entry)
            idx = int(np.argmin(distances))
            # return that entry's values directly (shape (n_value_fields,))
            return self._values[idx].copy()

        # pick k nearest neighbors (at most 4, or fewer if table is small)
        k = min(4, len(distances))
        # sort distances and take the k smallest indices
        nearest_idx = np.argsort(distances)[:k]
        # compute weights as 1/distance (closer entries get more influence)
        # 1e-10 prevents division by zero
        weights = 1.0 / (distances[nearest_idx] + 1e-10)
        # normalize weights so they sum to 1
        weights = weights/weights.sum()

        # weighted average of the k nearest neighbor values
        n_values = self._values.shape[1]
        interpolated_values = np.zeros(n_values)
        for i in range(k):
            interpolated_values += weights[i] * self._values[nearest_idx[i]]

        return interpolated_values


