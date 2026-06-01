"""SeconDose Core engine router.

Provides a factory function ``get_engine()`` that returns the appropriate
``DoseEngineBase`` implementation given a string key.  Engines are imported
lazily so unused engine modules (and their dependencies) are not loaded.

Engine keys
-----------
``"analytical"``
    Phase 1 ``SimpleAnalyticalDoseEngine`` wrapped in a thin adapter that
    satisfies the ``DoseEngineBase`` interface.  For regression testing only.
``"ccc"``
    Phase 2 ``CollapsedConeDoseEngine``.  Primary 2ndCheck production engine.

Reserved (not yet implemented):
    ``"monte_carlo"`` — Optional Phase 3+ GPU service.  Do not activate.

Usage
-----
    from DoseCalc.dose_engine.engine_router import get_engine

    engine = get_engine("ccc", kernel_path="kernels/6mv_placeholder.npz")
    dose   = engine.compute_dose(plan, ct, calibration)
"""
from __future__ import annotations

import importlib
from typing import Any

from .base import DoseEngineBase, DoseEngineError

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
# Values are (module_path, class_name) tuples for lazy import.
# Monte Carlo slot is reserved but commented out — do not uncomment in Phase 2.
_ENGINE_REGISTRY: dict[str, tuple[str, str]] = {
    "analytical": (
        "DoseCalc.dose_engine.analytical_adapter",
        "AnalyticalEngineAdapter",
    ),
    "ccc": (
        "DoseCalc.dose_engine.ccc_engine",
        "CollapsedConeDoseEngine",
    ),
    # "monte_carlo": (
    #     "DoseCalc.dose_engine.mc_engine",   # Phase 3+ only — not implemented
    #     "MonteCarloEngine",
    # ),
}

# Public list of valid keys for error messages.
VALID_ENGINE_KEYS: tuple[str, ...] = tuple(_ENGINE_REGISTRY.keys())


def get_engine(engine_key: str, **kwargs: Any) -> DoseEngineBase:
    """Instantiate and return the engine identified by *engine_key*.

    Parameters
    ----------
    engine_key:
        One of the registered keys (``"analytical"``, ``"ccc"``).
    **kwargs:
        Passed verbatim to the engine's ``__init__``.

    Returns
    -------
    DoseEngineBase
        A fully constructed engine instance ready for ``compute_dose``.

    Raises
    ------
    DoseEngineError
        If *engine_key* is not registered or instantiation fails.
    """
    key = str(engine_key).strip().lower()
    if key not in _ENGINE_REGISTRY:
        raise DoseEngineError(
            f"Unknown engine key {engine_key!r}. "
            f"Valid keys: {list(VALID_ENGINE_KEYS)}"
        )
    module_path, class_name = _ENGINE_REGISTRY[key]
    try:
        module = importlib.import_module(module_path)
    except ImportError as exc:
        raise DoseEngineError(
            f"Engine module '{module_path}' could not be imported: {exc}"
        ) from exc
    try:
        engine_cls = getattr(module, class_name)
    except AttributeError as exc:
        raise DoseEngineError(
            f"Engine class '{class_name}' not found in module '{module_path}': {exc}"
        ) from exc
    try:
        instance = engine_cls(**kwargs)
    except Exception as exc:
        raise DoseEngineError(
            f"Engine '{engine_key}' instantiation failed: {exc}"
        ) from exc

    if not isinstance(instance, DoseEngineBase):
        raise DoseEngineError(
            f"Engine '{engine_key}' ({class_name}) does not implement DoseEngineBase."
        )
    return instance

