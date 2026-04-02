from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import h5py
import numpy as np


@dataclass(frozen=True, slots=True)
class PassRecord:
    """A single calibrated satellite pass."""

    source_file: Path
    tile: str
    ref: str
    pointing: int
    pass_id: str
    norad_id: int
    unix_time: float
    noise_db: float
    alt_deg: np.ndarray
    az_deg: np.ndarray
    power_db: np.ndarray

    def __repr__(self) -> str:
        return (
            f"PassRecord(norad_id={self.norad_id}, "
            f"pointing={self.pointing}, "
            f"n_samples={len(self.power_db)}, "
            f"peak_alt={np.nanmax(self.alt_deg):.2f}°)"
        )


class PassFile:
    """Reader for a single EMBERS pass HDF5 file."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        if not self.path.is_file():
            raise FileNotFoundError(self.path)

        with h5py.File(self.path, "r") as h5:
            try:
                self.tile = str(h5.attrs["tile"])
                self.ref = str(h5.attrs["ref"])
            except KeyError as exc:
                raise ValueError(
                    f"{self.path} is missing required attribute: {exc}"
                ) from exc

            if "pointings" not in h5:
                raise ValueError(
                    f"{self.path} does not follow expected layout "
                    "(/pointings/<point>/passes/...)"
                )

            self.pointings = tuple(sorted(int(k) for k in h5["pointings"].keys()))

    def __repr__(self) -> str:
        return (
            f"PassFile(path={self.path!s}, tile={self.tile!r}, "
            f"ref={self.ref!r}, pointings={self.pointings})"
        )

    def _make_record(
        self,
        grp: h5py.Group,
        *,
        pointing: int,
        pass_id: str,
    ) -> PassRecord:
        """Construct a PassRecord from a single pass group."""
        return PassRecord(
            source_file=self.path,
            tile=self.tile,
            ref=self.ref,
            pointing=pointing,
            pass_id=pass_id,
            norad_id=int(grp.attrs["norad_id"]),
            unix_time=float(grp.attrs["unix_time"]),
            noise_db=float(grp.attrs["noise_db"]),
            alt_deg=np.asarray(grp["alt_deg"]),
            az_deg=np.asarray(grp["az_deg"]),
            power_db=np.asarray(grp["power_db"]),
        )

    def iter_passes(self, pointing: int) -> Iterator[PassRecord]:
        """Yield passes for a specific pointing."""
        if pointing not in self.pointings:
            raise ValueError(
                f"Pointing {pointing} not in {self.path.name}. "
                f"Available: {self.pointings}"
            )

        with h5py.File(self.path, "r") as h5:
            passes_group = h5["pointings"][str(pointing)]["passes"]

            for pass_id in sorted(passes_group.keys()):
                yield self._make_record(
                    passes_group[pass_id],
                    pointing=pointing,
                    pass_id=pass_id,
                )

    def read_passes(self, pointing: int) -> list[PassRecord]:
        """Return passes for a specific pointing."""
        return list(self.iter_passes(pointing))
