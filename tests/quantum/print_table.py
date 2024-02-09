import pandas as pd
import os, sys
from functools import reduce
from math import inf

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"

pd.set_option('display.max_rows', 500)

def process_averages(table: str, observable = "all", kappa = "all", skip_kappa = None, remove_outliers: bool = True, without_infinity: bool = False, add_times: bool = False, cpp: bool = False):
    data = pd.read_csv(os.path.join(SCRIPT_DIR, "results", f"[result{'-cpp' if cpp else ''}]{table}.csv"))

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
            
    if without_infinity and "time_lumping" in data.columns:
        data = pd.DataFrame(
            [row for (_,row) in data.iterrows() if (
                    row["time_lumping"] != inf
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

    return data.groupby(by=grouping).mean(numeric_only=True)

def table1(*, cpp:bool=False):
    if not cpp:
        time_grover_clue = process_averages("q_search_full_clue", skip_kappa=1, remove_outliers=False, add_times=True)["total time"].droplevel(-1).to_frame().rename(columns={"total time" : ("Grover", "CLUE")})
        time_grover_ddsim = process_averages("q_search_full_ddsim", skip_kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("Grover", "DDSIM")})
        time_sat_ddsim = process_averages("q_sat_full_ddsim", skip_kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("SAT", "DDSIM")})
        time_sat_clue = process_averages("q_sat_full_direct", skip_kappa=1, remove_outliers=False, add_times=True)["total time"].droplevel(-1).to_frame().rename(columns={"total time" : ("SAT", "CLUE")})
        d_sat_clue = process_averages("q_sat_direct", remove_outliers=False)["red. size"].to_frame().rename(columns={"red. size" : ("SAT", "d")})
        time_maxcut_ddsim = process_averages("q_maxcut_full_ddsim", skip_kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("MaxCut", "DDSIM")})
        time_maxcut_clue = process_averages("q_maxcut_full_direct", skip_kappa=1, remove_outliers=False, add_times=True)["total time"].droplevel(-1).to_frame().rename(columns={"total time" : ("MaxCut", "CLUE")})
        d_maxcut_clue = process_averages("q_maxcut_direct", remove_outliers=False)["red. size"].to_frame().rename(columns={"red. size" : ("MaxCut", "d")})

        return reduce(
            lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
            [time_grover_ddsim, time_grover_clue, time_sat_ddsim, time_sat_clue, d_sat_clue, time_maxcut_ddsim, time_maxcut_clue, d_maxcut_clue])
    else:
        time_grover_clue = process_averages("q_search_clue", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("Grover", "CLUE")})
        time_grover_ddsim = process_averages("q_search_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("Grover", "DDSIM")})
        time_grover_ddsim_alone = process_averages("q_search_full_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("Grover", "DDSIM-ALONE")})
        time_sat_clue = process_averages("q_sat_clue", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("SAT", "CLUE")})
        time_sat_ddsim = process_averages("q_sat_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("SAT", "DDSIM")})
        time_sat_ddsim_alone = process_averages("q_sat_full_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("SAT", "DDSIM-ALONE")})
        d_sat_clue = process_averages("q_sat_clue", remove_outliers=False, add_times=False, cpp=True)["red. ratio"].droplevel(0).to_frame().rename(columns={"red. ratio" : ("SAT", "d-CLUE")})
        for i in range(len(d_sat_clue)):
            d_sat_clue.iloc[i][0] = 2**d_sat_clue.iloc[i].name * d_sat_clue.iloc[i][0]
        d_sat_ddsim = process_averages("q_sat_ddsim", remove_outliers=False, add_times=False, cpp=True)["red. ratio"].droplevel(0).to_frame().rename(columns={"red. ratio" : ("SAT", "d-DD")})
        for i in range(len(d_sat_ddsim)):
            d_sat_ddsim.iloc[i][0] = 2**d_sat_ddsim.iloc[i].name * d_sat_ddsim.iloc[i][0]
        time_maxcut_clue = process_averages("q_maxcut_clue", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("MAXCUT", "CLUE")})
        time_maxcut_ddsim = process_averages("q_maxcut_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("MAXCUT", "DDSIM")})
        time_maxcut_ddsim_alone = process_averages("q_maxcut_full_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("MAXCUT", "DDSIM-ALONE")})
        d_maxcut_clue = process_averages("q_maxcut_clue", remove_outliers=False, add_times=False, cpp=True)["red. ratio"].droplevel(0).to_frame().rename(columns={"red. ratio" : ("MAXCUT", "d")})
        for i in range(len(d_maxcut_clue)):
            d_maxcut_clue.iloc[i][0] = 2**d_maxcut_clue.iloc[i].name * d_maxcut_clue.iloc[i][0]
        d_maxcut_ddsim = process_averages("q_maxcut_ddsim", remove_outliers=False, add_times=False, cpp=True)["red. ratio"].droplevel(0).to_frame().rename(columns={"red. ratio" : ("MAXCUT", "d-DD")})
        for i in range(len(d_maxcut_ddsim)):
            d_maxcut_ddsim.iloc[i][0] = 2**d_maxcut_ddsim.iloc[i].name * d_maxcut_ddsim.iloc[i][0]

        columns = [
                # Grover
                time_grover_clue, 
                time_grover_ddsim, 
                time_grover_ddsim_alone, 
                # SAT
                time_sat_clue, 
                time_sat_ddsim, 
                time_sat_ddsim_alone, 
                d_sat_clue, 
                d_sat_ddsim,
                # CUT
                time_maxcut_clue, 
                time_maxcut_ddsim, 
                time_maxcut_ddsim_alone, 
                d_maxcut_clue, 
                d_maxcut_ddsim
            ]
        
        # return columns

        return reduce(
            lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
            columns)


def table2(*, cpp:bool=False):
    red_zero = process_averages("q_benchmark_clue", remove_outliers=False, observable=0)["red. ratio"].droplevel(-1).to_frame().rename(columns={"red. ratio": "d/N wrt S_0"})
    avg_clue = process_averages("q_benchmark_clue", remove_outliers=False)[["red. ratio","time_lumping"]].rename(columns={"red. ratio": "Avg. d/N across S_x", "time_lumping": "Avg. time (s)"})
    avg_time_ddsim = process_averages("q_benchmark_full_ddsim", remove_outliers=False, kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration": "DDSIM time"})

    return reduce(
        lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
        [red_zero, avg_clue, avg_time_ddsim])

def satcpp():
    return process_averages("q_sat_clue", remove_outliers=False, add_times=False,cpp=True)

if __name__ == "__main__":
    what = sys.argv[1]
    if what == "table1":
        print(table1())
    elif what == "table2":
        print(table2())
    elif what == "table1-cpp":
        print(table1(cpp=True))
    elif what == "table2-cpp":
        print(table2(cpp=False))
    else:
        n = 2; observable = "all"; kappa = "all"; skip_kappa = None; remove_outliers = True; without_infinity = False; add_times = False; is_cpp = False
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
                elif sys.argv[n].endswith("cpp"):
                    is_cpp = True; n+=1
                else:
                    n += 1
            else:
                n += 1

        print(
            process_averages(what, observable=observable, kappa=kappa, skip_kappa=skip_kappa, remove_outliers=remove_outliers, without_infinity=without_infinity, add_times=add_times,cpp=is_cpp)
        )
