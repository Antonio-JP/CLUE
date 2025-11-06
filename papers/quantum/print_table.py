import pandas as pd
import os, sys
from functools import reduce
from math import inf
import numpy as np

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"

pd.set_option('display.max_rows', 500)

def process_averages(table: str, observable = "all", kappa = "all", skip_kappa = None, remove_outliers: bool = True, without_infinity: bool = False, add_times: bool = False):
    data = pd.read_csv(os.path.join(SCRIPT_DIR, "results", f"[result]{table}.csv"))

    ## CLEANING THE DATA
    if remove_outliers and len(data) > 20: # we require enough data to remove outliers
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
                )], columns=data.columns)
            
    if "red. ratio" in data.columns:
        data.insert(list(data.columns).index("red. ratio") + 1, "red. size", pd.Series([2**row["size"] * (float(row["red. ratio"]) if row["red. ratio"] != "unknown" else inf) for (_,row) in data.iterrows()]))
        data["red. ratio"] = pd.Series([(float(row["red. ratio"]) if row["red. ratio"] != "unknown" else inf) for (_,row) in data.iterrows()])
    
    if add_times and len([col for col in data.columns if "time" in col]) > 1:
        data.insert(len(data.columns), "total time", pd.Series([sum(float(row[col]) for col in data.columns if col.startswith("time")) for (_,row) in data.iterrows()]))

    ## FILTERING BY OBSERVABLE IF REQUIRED
    grouping = (["name"] if "name" in data.columns else []) + ["size"] + (["obs"] if ("obs" in data.columns and observable != "all") else []) + (["kappa"] if "kappa" in data.columns else [])

    if "obs" in data.columns:
        if observable == "all": # we simply remove the column
            data = data.groupby(grouping + ["obs"]).mean(numeric_only=True)
            data = data.reset_index()
            data = data.drop("obs", axis=1)
        elif observable != "split":
            data = data[data["obs"] == str(observable)]

    if "kappa" in data.columns:
        if kappa != "all":
            data = data[(data["kappa"] == kappa) & (data["kappa"] != skip_kappa)]
        else:
            data = data[data["kappa"] != skip_kappa]
    
    inf_data = None
    if without_infinity and any(col in data.columns for col in ["time_lumping", "time_iteration"]):
        col = "time_lumping" if "time_lumping" in data.columns else "time_iteration"
        inf_data = pd.DataFrame([row for (_,row) in data.iterrows() if row[col] == inf], columns=data.columns)
        data = pd.DataFrame([row for (_,row) in data.iterrows() if row[col] != inf], columns=data.columns)

    result = data.groupby(by=grouping).mean(numeric_only=True)
    
    if inf_data is not None:
        inf_data = inf_data.groupby(by=grouping).count()[col].to_frame().rename(columns={col: "inf. count"})
        result = result.join(inf_data, how="outer")
    return result

def table1():
    grover_clue = process_averages("q_search_full_clue", kappa=1, remove_outliers=False, add_times=True, without_infinity=True)
    time_grover_clue = grover_clue["total time"].droplevel(-1).to_frame().rename(columns={"total time" : ("Grover", "CLUE")})
    inf_time_grover_clue = grover_clue["inf. count"].droplevel(-1).to_frame().rename(columns={"inf. count" : ("Grover", "CL-∞")})
    
    grover_ddsim = process_averages("q_search_full_ddsim", skip_kappa=1, without_infinity=True)
    time_grover_ddsim = grover_ddsim["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("Grover", "DDSIM")})
    inf_time_grover_ddsim = grover_ddsim["inf. count"].droplevel(-1).to_frame().rename(columns={"inf. count" : ("Grover", "DD-∞")})
    
    grover_quokka = process_averages("q_search_full_quokka#", skip_kappa=1, without_infinity=True)
    if not "time_iteration" in grover_quokka.columns:
        grover_quokka["time_iteration"] = np.nan
    time_grover_quokka = grover_quokka["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("Grover", "QUOKKA")})
    inf_time_grover_quokka = grover_quokka["inf. count"].droplevel(-1).to_frame().rename(columns={"inf. count" : ("Grover", "QK-∞")})
    
    sat_ddsim = process_averages("q_sat_full_ddsim", skip_kappa=1, without_infinity=True)
    time_sat_ddsim = sat_ddsim["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("SAT", "DDSIM")})
    inf_time_sat_ddsim = sat_ddsim["inf. count"].droplevel(-1).to_frame().rename(columns={"inf. count" : ("SAT", "DD-∞")})
    
    sat_quokka = process_averages("q_sat_full_quokka#", skip_kappa=1, without_infinity=True)
    time_sat_quokka = sat_quokka["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("SAT", "QUOKKA")})
    inf_time_sat_quokka = sat_quokka["inf. count"].droplevel(-1).to_frame().rename(columns={"inf. count" : ("SAT", "QK-∞")})
    
    sat_clue = process_averages("q_sat_full_direct", skip_kappa=1, remove_outliers=False, add_times=True, without_infinity=True)
    time_sat_clue = sat_clue["total time"].droplevel(-1).to_frame().rename(columns={"total time" : ("SAT", "CLUE")})
    inf_time_sat_clue = sat_clue["inf. count"].droplevel(-1).to_frame().rename(columns={"inf. count" : ("SAT", "CL-∞")})
    
    d_sat_clue = process_averages("q_sat_direct", remove_outliers=False)["red. size"].to_frame().rename(columns={"red. size" : ("SAT", "d")})
    
    maxcut_ddsim = process_averages("q_maxcut_full_ddsim", skip_kappa=1, without_infinity=True)
    time_maxcut_ddsim = maxcut_ddsim["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("MaxCut", "DDSIM")})
    inf_time_maxcut_ddsim = maxcut_ddsim["inf. count"].droplevel(-1).to_frame().rename(columns={"inf. count" : ("MaxCut", "DD-∞")})
    
    maxcut_quokka = process_averages("q_maxcut_full_quokka#", skip_kappa=1, without_infinity=True)
    time_maxcut_quokka = maxcut_quokka["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("MaxCut", "QUOKKA")})
    inf_time_maxcut_quokka = maxcut_quokka["inf. count"].droplevel(-1).to_frame().rename(columns={"inf. count" : ("MaxCut", "QK-∞")})
    
    maxcut_clue = process_averages("q_maxcut_full_direct", skip_kappa=1,remove_outliers=False, add_times=True, without_infinity=True)
    time_maxcut_clue = maxcut_clue["total time"].droplevel(-1).to_frame().rename(columns={"total time" : ("MaxCut", "CLUE")})
    inf_time_maxcut_clue = maxcut_clue["inf. count"].droplevel(-1).to_frame().rename(columns={"inf. count" : ("MaxCut", "CL-∞")})
    
    d_maxcut_clue = process_averages("q_maxcut_direct", remove_outliers=False)["red. size"].to_frame().rename(columns={"red. size" : ("MaxCut", "d")})

    return reduce(
        lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
        [time_grover_ddsim, inf_time_grover_ddsim, time_grover_quokka, inf_time_grover_quokka, time_grover_clue, inf_time_grover_clue, 
         time_sat_ddsim, inf_time_sat_ddsim, time_sat_quokka, inf_time_sat_quokka, time_sat_clue, inf_time_sat_clue, d_sat_clue, 
         time_maxcut_ddsim, inf_time_maxcut_ddsim, time_maxcut_quokka, inf_time_maxcut_quokka, time_maxcut_clue, inf_time_maxcut_clue, d_maxcut_clue])

def table2():
    red_zero = process_averages("q_benchmark_clue", remove_outliers=False, observable=0)["red. ratio"].droplevel(-1).to_frame().rename(columns={"red. ratio": "d/N wrt S_0"})
    avg_clue = process_averages("q_benchmark_clue", remove_outliers=False)[["red. ratio","time_lumping"]].rename(columns={"red. ratio": "Avg. d/N across S_x", "time_lumping": "Avg. time (s)"})
    avg_time_ddsim = process_averages("q_benchmark_full_ddsim", remove_outliers=False, kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration": "DDSIM time"})
    avg_time_quokka = process_averages("q_benchmark_full_quokka#", remove_outliers=False, kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration": "QUOKKA time"})

    return reduce(
        lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
        [red_zero, avg_clue, avg_time_ddsim, avg_time_quokka])

if __name__ == "__main__":
    what = sys.argv[1]
    if what == "table1":
        print(table1())
    elif what == "table2":
        print(table2())
    else:
        n = 2; observable = "all"; kappa = "all"; skip_kappa = None; remove_outliers = True; without_infinity = False; add_times = False
        while n < len(sys.argv):
            if sys.argv[n].startswith("-"):
                if sys.argv[n].endswith("obs"):
                    observable = sys.argv[n+1]; n += 2
                elif sys.argv[n].endswith("nk"):
                    skip_kappa = int(sys.argv[n+1]); n += 2
                elif sys.argv[n].endswith("k"):
                    kappa = int(sys.argv[n+1]); n += 2
                elif sys.argv[n].endswith("wo"):
                    remove_outliers = False; n+=1
                elif sys.argv[n].endswith("noinf"):
                    without_infinity = True; n+=1
                elif sys.argv[n].endswith("add"):
                    add_times = True; n+=1
                else:
                    n += 1
            else:
                n += 1

        print(process_averages(what, observable=observable, kappa=kappa, skip_kappa=skip_kappa, remove_outliers=remove_outliers, without_infinity=without_infinity, add_times=add_times))
