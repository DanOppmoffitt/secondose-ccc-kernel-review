from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[1]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from DoseCalc.scripts.fit_experimental_field_size_hybrid_kernel import main


if __name__ == "__main__":
    raise SystemExit(main())

