"""Download the open GPS sample data used by the experiments and tests.

Station 0759 (GSI, Japan), 2005-04-02, 1 hour of dual-frequency GPS at 30 s,
with the matching broadcast navigation file, as distributed with RTKLIB.
"""

import urllib.request
from pathlib import Path

BASE = "https://raw.githubusercontent.com/tomojitakasu/RTKLIB/master/test/data/rinex/"
FILES = ["07590920.05o", "07590920.05n"]
DATA = Path(__file__).resolve().parents[1] / "data"


def main():
    DATA.mkdir(exist_ok=True)
    for name in FILES:
        target = DATA / name
        if target.exists():
            print(f"{name}: already present")
            continue
        urllib.request.urlretrieve(BASE + name, target)
        print(f"{name}: downloaded ({target.stat().st_size / 1024:.0f} kB)")


if __name__ == "__main__":
    main()
