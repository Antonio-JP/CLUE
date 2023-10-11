import pandas as pd
import os, sys

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"

if __name__ == "__main__":
    what = sys.argv[1]
    data = pd.read_csv(os.path.join(SCRIPT_DIR, f"[result]{what}.csv"))
    ## CLEANING THE DATA
    if "time_lumping" in data.columns: ## We remove outlier data
        grouped = data.groupby(by=["size"] + (["kappa"] if "kappa" in data.columns else []))
        q_low = grouped["time_lumping"].quantile(0.05)
        q_hi  = grouped["time_lumping"].quantile(0.95)
        M_low = grouped["memory (MB)"].quantile(0.05)
        M_high = grouped["memory (MB)"].quantile(0.95)
        get_filter = (lambda D, r : D[r["size"]]) if (not "kappa" in data.columns) else (lambda D, r : D[r["size"],r["kappa"]])
        data = pd.DataFrame(
            [row for (_,row) in data.iterrows() if (
                row["time_lumping"] < get_filter(q_hi, row) and 
                row["time_lumping"] > get_filter(q_low,row) and
                row["memory (MB)"] < get_filter(M_high, row) and 
                row["memory (MB)"] > get_filter(M_low, row)
            )])
        
        if "red. ratio" in data.columns:
            data["red. size"] = (2**data["size"])*data["red. ratio"]
    
    ## PRINTING RESULTING DATA
    print(data.drop(columns=data.columns[-2:-1]).groupby(by=["size"] + (["kappa"] if "kappa" in data.columns else [])).mean())
    sys.exit(1)