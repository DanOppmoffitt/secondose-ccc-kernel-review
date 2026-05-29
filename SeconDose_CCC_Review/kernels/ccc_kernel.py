"""CCC energy deposition kernel data structures and I/O.

A ``CCCKernelData`` object represents a pre-computed 3-D energy deposition
kernel (EDK) in polar form: ``K[energy_bin, r_index, theta_index]`` or
``K[r_index, theta_index]`` for a single effective energy.

Storage format: NumPy ``.npz`` compressed archive with a JSON metadata
sidecar embedded under the key ``"_meta_json"``.

Integrity: The kernel matrix bytes are SHA-256 hashed at save time.  The
hash is stored in ``_meta_json`` and re-verified at load time.  A checksum
mismatch raises ``CCCKernelIntegrityError``.

Placeholder kernel: ``build_placeholder_ccc_kernel`` returns a physically
plausible (but not validated) kernel suitable for infrastructure testing.
The placeholder uses a single effective energy bin with an exponential-decay
primary component and a Gaussian scatter tail, scaled to integrate to a
realistic deposited fraction (~0.95 for a large water phantom).
"""
from __future__ import annotations

import hashlib
import io
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

_log = logging.getLogger(__name__)

# Reasonable bounds for the deposited-energy fraction.
# Total absorbed dose / total TERMA for a 10x10 cm field in a large water
# phantom is typically 0.92-0.98.  Values outside [0.5, 1.0] indicate a
# malformed kernel.
_DEPOSITED_FRACTION_BOUNDS = (0.50, 1.00)

# Maximum radial extent of the kernel (cm).  Kernels with r_grid_cm extending
# beyond this value are flagged as suspicious.
_MAX_KERNEL_RADIUS_CM = 35.0


class CCCKernelError(Exception):
    """Base for CCC kernel errors."""


class CCCKernelIntegrityError(CCCKernelError):
    """Raised when SHA-256 checksum of loaded kernel does not match stored value."""


@dataclass(frozen=True)
class CCCKernelData:
    """Immutable polar energy deposition kernel for collapsed-cone convolution.

    Attributes
    ----------
    source_citation : str
        Bibliographic reference or source ID (e.g. ``"Mackie1988"``).
    energy_bins_mev : np.ndarray, shape (N_bins,)
        Photon energies in MeV for each spectral component.  A single-bin
        array represents a monoenergetic effective-energy model.
    fluence_weights : np.ndarray, shape (N_bins,)
        Relative fluence weight for each energy bin; must sum to 1.0.
    r_grid_cm : np.ndarray, shape (N_r,)
        Radial distance grid in cm from the interaction point.  Must start at
        0 and be monotonically increasing.
    theta_grid_deg : np.ndarray, shape (N_theta,)
        Polar angle grid in degrees (0 = forward along beam direction,
        180 = backward).
    kernel_matrix : np.ndarray
        Energy deposition values.  Shape is either ``(N_r, N_theta)`` for a
        monoenergetic kernel or ``(N_bins, N_r, N_theta)`` for polyenergetic.
        Units: fraction of TERMA deposited per (cm·sr) bin.  The full
        integral over all r and theta equals ``deposited_fraction``.
    deposited_fraction : float
        Expected value of integral(K) for a large water phantom.  Typically
        0.92–0.98 for a 6 MV beam.  Used as a sanity-check sentinel.
    created_date : str
        ISO-8601 date of kernel creation (``"YYYY-MM-DD"``).
    checksum : str
        SHA-256 hex digest of the raw ``kernel_matrix`` bytes
        (``kernel_matrix.tobytes()``).  Verified at load time.
    notes : str
        Optional free-text provenance notes.
    """

    source_citation: str
    energy_bins_mev: np.ndarray
    fluence_weights: np.ndarray
    r_grid_cm: np.ndarray
    theta_grid_deg: np.ndarray
    kernel_matrix: np.ndarray
    deposited_fraction: float
    created_date: str
    checksum: str
    notes: str = ""

    def __post_init__(self) -> None:
        e = np.asarray(self.energy_bins_mev, dtype=np.float64)
        w = np.asarray(self.fluence_weights, dtype=np.float64)
        r = np.asarray(self.r_grid_cm, dtype=np.float64)
        theta = np.asarray(self.theta_grid_deg, dtype=np.float64)
        k = np.asarray(self.kernel_matrix, dtype=np.float64)

        if e.ndim != 1 or e.size == 0:
            raise CCCKernelError("energy_bins_mev must be a non-empty 1-D array.")
        if w.ndim != 1 or w.size != e.size:
            raise CCCKernelError("fluence_weights must match energy_bins_mev in length.")
        if not np.isclose(float(np.sum(w)), 1.0, atol=1e-6):
            raise CCCKernelError(
                f"fluence_weights must sum to 1.0; got {float(np.sum(w)):.8f}."
            )
        if r.ndim != 1 or r.size < 2:
            raise CCCKernelError("r_grid_cm must be a 1-D array with at least 2 points.")
        if not np.all(np.diff(r) > 0):
            raise CCCKernelError("r_grid_cm must be monotonically increasing.")
        if float(r[0]) < 0.0:
            raise CCCKernelError("r_grid_cm must start at >= 0.")
        if theta.ndim != 1 or theta.size < 2:
            raise CCCKernelError("theta_grid_deg must be a 1-D array with at least 2 points.")

        n_bins = int(e.size)
        n_r = int(r.size)
        n_theta = int(theta.size)
        if n_bins == 1:
            expected_shape = (n_r, n_theta)
        else:
            expected_shape = (n_bins, n_r, n_theta)
        if k.shape != expected_shape:
            raise CCCKernelError(
                f"kernel_matrix shape {k.shape} does not match expected {expected_shape}."
            )
        if not np.all(np.isfinite(k)):
            raise CCCKernelError("kernel_matrix contains non-finite values.")
        if np.any(k < 0.0):
            raise CCCKernelError("kernel_matrix contains negative values.")

        df = float(self.deposited_fraction)
        if not (_DEPOSITED_FRACTION_BOUNDS[0] <= df <= _DEPOSITED_FRACTION_BOUNDS[1]):
            raise CCCKernelError(
                f"deposited_fraction {df:.4f} is outside reasonable bounds "
                f"{_DEPOSITED_FRACTION_BOUNDS}."
            )

        # Freeze arrays
        object.__setattr__(self, "energy_bins_mev", e)
        object.__setattr__(self, "fluence_weights", w)
        object.__setattr__(self, "r_grid_cm", r)
        object.__setattr__(self, "theta_grid_deg", theta)
        object.__setattr__(self, "kernel_matrix", k)
        object.__setattr__(self, "deposited_fraction", df)

    @property
    def n_bins(self) -> int:
        return int(self.energy_bins_mev.size)

    @property
    def n_r(self) -> int:
        return int(self.r_grid_cm.size)

    @property
    def n_theta(self) -> int:
        return int(self.theta_grid_deg.size)

    @property
    def is_monoenergetic(self) -> bool:
        return self.n_bins == 1


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_ccc_kernel(kernel: CCCKernelData) -> list[str]:
    """Return a list of warning / issue strings for *kernel*.

    An empty list means no issues found.  This function does not raise;
    callers can decide how to handle warnings.
    """
    issues: list[str] = []

    if float(kernel.r_grid_cm[-1]) > _MAX_KERNEL_RADIUS_CM:
        issues.append(
            f"r_grid_cm extends to {float(kernel.r_grid_cm[-1]):.1f} cm, "
            f"which exceeds the expected maximum of {_MAX_KERNEL_RADIUS_CM} cm."
        )
    if float(kernel.r_grid_cm[0]) != 0.0:
        issues.append(
            f"r_grid_cm[0] = {float(kernel.r_grid_cm[0]):.4f} cm; expected 0.0."
        )
    if kernel.is_monoenergetic:
        issues.append(
            "Kernel is monoenergetic (single energy bin).  "
            "Polyenergetic representation is recommended for 6 MV (Decision A)."
        )
    if float(np.sum(kernel.kernel_matrix)) == 0.0:
        issues.append("kernel_matrix integral is zero — likely a malformed kernel.")

    # Verify checksum
    computed = _compute_checksum(kernel.kernel_matrix)
    if computed != kernel.checksum:
        issues.append(
            f"Checksum mismatch: stored={kernel.checksum}, computed={computed}."
        )

    return issues


# ---------------------------------------------------------------------------
# Save / Load
# ---------------------------------------------------------------------------

def _compute_checksum(matrix: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(matrix, dtype=np.float64).tobytes()).hexdigest()


def save_ccc_kernel(kernel: CCCKernelData, path: str | Path) -> None:
    """Save *kernel* to a ``.npz`` file with embedded JSON metadata.

    Parameters
    ----------
    kernel:
        A validated ``CCCKernelData`` instance.
    path:
        Destination file path (``.npz`` extension recommended).
    """
    p = Path(path)
    meta: dict[str, Any] = {
        "source_citation": kernel.source_citation,
        "deposited_fraction": float(kernel.deposited_fraction),
        "created_date": kernel.created_date,
        "checksum": kernel.checksum,
        "notes": kernel.notes,
        "n_bins": kernel.n_bins,
        "n_r": kernel.n_r,
        "n_theta": kernel.n_theta,
        "is_monoenergetic": kernel.is_monoenergetic,
    }
    meta_bytes = json.dumps(meta, indent=2).encode("utf-8")
    np.savez_compressed(
        p,
        _meta_json=np.frombuffer(meta_bytes, dtype=np.uint8),
        energy_bins_mev=kernel.energy_bins_mev,
        fluence_weights=kernel.fluence_weights,
        r_grid_cm=kernel.r_grid_cm,
        theta_grid_deg=kernel.theta_grid_deg,
        kernel_matrix=kernel.kernel_matrix,
    )
    _log.info("Saved CCC kernel to %s (checksum=%s)", p, kernel.checksum)


def load_ccc_kernel(path: str | Path) -> CCCKernelData:
    """Load and verify a CCC kernel from *path*.

    Raises
    ------
    CCCKernelIntegrityError
        If the stored checksum does not match the loaded kernel matrix.
    CCCKernelError
        If required arrays are missing from the file.
    """
    p = Path(path)
    data = np.load(p, allow_pickle=False)

    required = {"_meta_json", "energy_bins_mev", "fluence_weights",
                "r_grid_cm", "theta_grid_deg", "kernel_matrix"}
    missing = required - set(data.files)
    if missing:
        raise CCCKernelError(f"Kernel file {p} is missing arrays: {missing}")

    meta: dict[str, Any] = json.loads(data["_meta_json"].tobytes().decode("utf-8"))

    kernel_matrix = np.asarray(data["kernel_matrix"], dtype=np.float64)
    computed_checksum = _compute_checksum(kernel_matrix)
    stored_checksum = str(meta.get("checksum", ""))
    if computed_checksum != stored_checksum:
        raise CCCKernelIntegrityError(
            f"Kernel integrity check failed for {p}.\n"
            f"  stored  checksum: {stored_checksum}\n"
            f"  computed checksum: {computed_checksum}\n"
            "The file may be corrupted or have been modified."
        )

    return CCCKernelData(
        source_citation=str(meta["source_citation"]),
        energy_bins_mev=np.asarray(data["energy_bins_mev"], dtype=np.float64),
        fluence_weights=np.asarray(data["fluence_weights"], dtype=np.float64),
        r_grid_cm=np.asarray(data["r_grid_cm"], dtype=np.float64),
        theta_grid_deg=np.asarray(data["theta_grid_deg"], dtype=np.float64),
        kernel_matrix=kernel_matrix,
        deposited_fraction=float(meta["deposited_fraction"]),
        created_date=str(meta["created_date"]),
        checksum=stored_checksum,
        notes=str(meta.get("notes", "")),
    )


# ---------------------------------------------------------------------------
# Placeholder kernel for infrastructure testing
# ---------------------------------------------------------------------------

def build_placeholder_ccc_kernel(
    *,
    n_r: int = 60,
    r_max_cm: float = 30.0,
    n_theta: int = 48,
    primary_decay_cm: float = 7.0,
    scatter_sigma_cm: float = 4.0,
    scatter_weight: float = 0.18,
    deposited_fraction: float = 0.95,
    energy_mev: float = 1.75,
    source_citation: str = "placeholder_6MV_infrastructure_test",
    notes: str = "Non-validated placeholder kernel for Phase 2 infrastructure tests.",
) -> CCCKernelData:
    """Build a physically plausible but NOT validated CCC kernel.

    The kernel uses a single effective energy bin representing the mean photon
    energy of a 6 MV beam (~1.75 MeV).  The radial shape is a two-component
    model: exponential primary decay + Gaussian scatter tail.

    This kernel MUST NOT be used for any clinical or research dose calculation.
    It exists solely to allow testing of kernel I/O, engine construction,
    and infrastructure code before real measured kernels are available.

    Parameters
    ----------
    n_r:
        Number of radial grid points.
    r_max_cm:
        Maximum radial distance in cm.
    n_theta:
        Number of angular bins (should match Decision G: 24, 48, or 96).
    primary_decay_cm:
        Exponential decay constant for the primary kernel component.
    scatter_sigma_cm:
        Gaussian sigma for the scatter kernel component.
    scatter_weight:
        Weight of the scatter component (primary weight = 1 - scatter_weight).
    deposited_fraction:
        Target absorbed fraction (see ``CCCKernelData.deposited_fraction``).
    energy_mev:
        Effective photon energy in MeV for the single spectral bin.
    source_citation, notes:
        Provenance strings.
    """
    r = np.linspace(0.0, r_max_cm, n_r)
    theta_deg = np.linspace(0.0, 180.0, n_theta)

    # Forward-peaked primary component (exponential in r, cosine in theta)
    R, THETA = np.meshgrid(r, theta_deg, indexing="ij")
    theta_rad = np.deg2rad(THETA)
    forward_weight = 0.5 * (1.0 + np.cos(theta_rad))  # 1 at 0°, 0 at 180°
    primary = (1.0 - scatter_weight) * forward_weight * np.exp(-R / primary_decay_cm)

    # Isotropic scatter component
    scatter = scatter_weight * np.exp(-0.5 * (R / scatter_sigma_cm) ** 2)

    raw = primary + scatter
    raw = np.maximum(raw, 0.0)

    # Normalise so that the kernel integral equals deposited_fraction.
    total = float(np.sum(raw))
    if total <= 0.0:
        raise CCCKernelError("Placeholder kernel integral is zero — check parameters.")
    scale = deposited_fraction / total
    kernel_matrix = (raw * scale).astype(np.float64)

    energy_bins = np.array([energy_mev], dtype=np.float64)
    fluence_weights = np.array([1.0], dtype=np.float64)
    checksum = _compute_checksum(kernel_matrix)

    return CCCKernelData(
        source_citation=source_citation,
        energy_bins_mev=energy_bins,
        fluence_weights=fluence_weights,
        r_grid_cm=r,
        theta_grid_deg=theta_deg,
        kernel_matrix=kernel_matrix,
        deposited_fraction=deposited_fraction,
        created_date="2026-05-23",
        checksum=checksum,
        notes=notes,
    )

