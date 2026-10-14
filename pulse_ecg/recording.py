"""
The one signal type. PTB-XL, LUDB, QTDB and the device all become a Recording,
and nothing downstream knows or cares which one it came from.
"""

from dataclasses import dataclass, field

import numpy as np

STANDARD_LEADS = ("I", "II", "III", "aVR", "aVL", "aVF",
                  "V1", "V2", "V3", "V4", "V5", "V6")
_CANONICAL = {name.upper(): name for name in STANDARD_LEADS}

SOURCES = ("ptbxl", "ludb", "qtdb", "device", "synthetic")


def canonical_lead(name: str) -> str:
    """'AVF', 'avf' -> 'aVF'. Non-standard names (QTDB's 'MLII') pass through."""
    name = str(name).strip()
    return _CANONICAL.get(name.upper(), name)


@dataclass(frozen=True)
class Recording:
    signal: np.ndarray        # (n_samples, n_leads), millivolts
    fs: int                   # samples per second
    leads: tuple[str, ...]    # ("I", ..., "V6")  or  ("II",) for the device
    duration_s: float
    source: str               # one of SOURCES
    record_id: str
    meta: dict = field(default_factory=dict, compare=False)  # age, sex, patient_id, ...

    def __post_init__(self):
        sig = np.asarray(self.signal, dtype=float)
        if sig.ndim != 2:
            raise ValueError(f"{self.record_id}: signal must be 2-D (n_samples, n_leads), "
                             f"got shape {sig.shape}")
        leads = tuple(canonical_lead(l) for l in self.leads)
        if len(leads) != sig.shape[1]:
            raise ValueError(f"{self.record_id}: {sig.shape[1]} signal columns "
                             f"but {len(leads)} lead names")
        if len(set(leads)) != len(leads):
            raise ValueError(f"{self.record_id}: duplicate lead names {leads}")
        if not (isinstance(self.fs, (int, np.integer)) and self.fs > 0):
            raise ValueError(f"{self.record_id}: fs must be a positive integer, got {self.fs!r}")
        if self.source not in SOURCES:
            raise ValueError(f"{self.record_id}: unknown source {self.source!r}")
        expected = sig.shape[0] / self.fs
        if abs(self.duration_s - expected) > 1.0 / self.fs:
            raise ValueError(f"{self.record_id}: duration_s={self.duration_s} but "
                             f"{sig.shape[0]} samples at {self.fs} Hz is {expected:.3f} s")
        object.__setattr__(self, "signal", sig)
        object.__setattr__(self, "leads", leads)
        object.__setattr__(self, "fs", int(self.fs))

    @classmethod
    def from_array(cls, signal, fs, leads, source, record_id, meta=None):
        """Build a Recording, deriving duration_s from the sample count."""
        sig = np.asarray(signal, dtype=float)
        return cls(sig, int(fs), tuple(leads), sig.shape[0] / fs, source,
                   str(record_id), dict(meta or {}))

    @property
    def n_samples(self) -> int:
        return self.signal.shape[0]

    @property
    def n_leads(self) -> int:
        return self.signal.shape[1]

    def has_leads(self, *names: str) -> bool:
        return all(canonical_lead(n) in self.leads for n in names)

    def lead(self, name: str) -> np.ndarray:
        name = canonical_lead(name)
        try:
            return self.signal[:, self.leads.index(name)]
        except ValueError:
            raise KeyError(f"{self.record_id}: lead {name} not present "
                           f"(have {', '.join(self.leads)})") from None

    def select(self, leads) -> "Recording":
        """A copy restricted to `leads` — e.g. ["II"] to emulate the device."""
        names = [canonical_lead(l) for l in leads]
        cols = [self.leads.index(n) for n in names] if self.has_leads(*names) else None
        if cols is None:
            missing = [n for n in names if n not in self.leads]
            raise KeyError(f"{self.record_id}: cannot select missing leads {missing}")
        return Recording(self.signal[:, cols].copy(), self.fs, tuple(names),
                         self.duration_s, self.source, self.record_id, dict(self.meta))
