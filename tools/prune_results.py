"""
Delete experiment result artifacts under Results/ to keep the repo clean.
Safe by default: refuses to delete if paths are outside Results/.
"""
from pathlib import Path
import shutil

RESULTS_ROOT = Path("Results")


def prune_results() -> None:
    if not RESULTS_ROOT.exists():
        print("Results/ does not exist; nothing to prune.")
        return
    for child in RESULTS_ROOT.iterdir():
        if child.is_dir():
            print(f"Removing directory {child}")
            shutil.rmtree(child, ignore_errors=True)
        elif child.is_file():
            print(f"Removing file {child}")
            try:
                child.unlink()
            except FileNotFoundError:
                pass


if __name__ == "__main__":
    prune_results()
