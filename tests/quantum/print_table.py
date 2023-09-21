import pandas as pd
import os, sys

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"

if __name__ == "__main__":
    what = sys.argv[1]
    data = pd.read_csv(os.path.join(SCRIPT_DIR, f"[result]{what}.csv"))
    print(data.drop(columns=data.columns[-1:]).groupby(by=["size"]).mean())
    sys.exit(1)