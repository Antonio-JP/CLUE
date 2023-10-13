import pandas as pd
import os, sys

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"

if __name__ == "__main__":
    what = sys.argv[1]
    data = pd.read_csv(os.path.join(SCRIPT_DIR, "results", f"[result]{what}.csv"))
    n = 2; observable = "all"; rem_outliers = True
    while n < len(sys.argv):
        if sys.argv[n].startswith("-"):
            if sys.argv[n].endswith("obs"):
                observable = sys.argv[n+1]; n += 2
            elif sys.argv[n].endswith("wo"):
                rem_outliers = False; n+=1
            else:
                n += 1
        else:
            n += 1
    ## CLEANING THE DATA
    if rem_outliers:
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

    ## FILTERING BY OBSERVABLE IF REQUIRED
    if "obs" in data.columns:
        if observable == "all": # we simply remove the column
            data = data.drop("obs", axis=1)
        elif observable != "split":
            data = pd.DataFrame([row for (_,row) in data.iterrows() if row["obs"] == observable], columns=data.columns)

    
    ## PRINTING RESULTING DATA
    print(data.drop(columns=data.columns[-2:-1]).groupby(by=(["name"] if "name" in data.columns else []) + ["size"] + (["obs"] if ("obs" in data.columns and observable != "all") else []) + (["kappa"] if "kappa" in data.columns else [])).mean(numeric_only=True))
    #sys.exit(1)
