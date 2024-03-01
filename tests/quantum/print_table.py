import pandas as pd
import os, sys
from functools import reduce
from math import inf

SCRIPT_DIR = os.path.dirname(__file__)

pd.set_option('display.max_rows', 500)

#######################################################################################
### PROCESSING METHODS
#######################################################################################
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
            if data.dtypes['obs'] == 'O': data = data[data["obs"] == str(observable)]
            else: data = data[data["obs"] == data.dtypes["obs"].type(observable)]

    if "kappa" in data.columns:
        if kappa != "all":
            data = data[(data["kappa"] == kappa) & (data["kappa"] != skip_kappa)]
        else:
            data = data[data["kappa"] != skip_kappa]

    return data.groupby(by=grouping).mean(numeric_only=True)

def table1(*, cpp:bool=False, out:str="text", formt="sci"):
    ## Creating the table (depends on C++ argument)
    if not cpp:
        time_grover_clue = process_averages("q_search_full_clue", skip_kappa=1, remove_outliers=False, add_times=True)["total time"].droplevel(-1).to_frame().rename(columns={"total time" : ("Grover", "CLUE")})
        time_grover_ddsim = process_averages("q_search_full_ddsim", skip_kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("Grover", "DDSIM+CLUE")})
        time_sat_ddsim = process_averages("q_sat_full_ddsim", skip_kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("SAT", "DDSIM+CLUE")})
        time_sat_clue = process_averages("q_sat_full_direct", skip_kappa=1, remove_outliers=False, add_times=True)["total time"].droplevel(-1).to_frame().rename(columns={"total time" : ("SAT", "CLUE")})
        d_sat_clue = process_averages("q_sat_direct", remove_outliers=False)["red. size"].to_frame().rename(columns={"red. size" : ("SAT", "d")})
        time_maxcut_ddsim = process_averages("q_maxcut_full_ddsim", skip_kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration" : ("MaxCut", "DDSIM+CLUE")})
        time_maxcut_clue = process_averages("q_maxcut_full_direct", skip_kappa=1, remove_outliers=False, add_times=True)["total time"].droplevel(-1).to_frame().rename(columns={"total time" : ("MaxCut", "CLUE")})
        d_maxcut_clue = process_averages("q_maxcut_direct", remove_outliers=False)["red. size"].to_frame().rename(columns={"red. size" : ("MaxCut", "d")})

        data = reduce(
            lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
            [time_grover_ddsim, time_grover_clue, time_sat_ddsim, time_sat_clue, d_sat_clue, time_maxcut_ddsim, time_maxcut_clue, d_maxcut_clue])
    else:
        time_grover_clue = process_averages("q_search_clue", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("Grover", "CLUE")})
        time_grover_ddsim = process_averages("q_search_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("Grover", "DDSIM+CLUE")})
        time_grover_ddsim_alone = process_averages("q_search_full_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("Grover", "DDSIM")})
        time_sat_clue = process_averages("q_sat_clue", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("SAT", "CLUE")})
        time_sat_ddsim = process_averages("q_sat_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("SAT", "DDSIM+CLUE")})
        time_sat_ddsim_alone = process_averages("q_sat_full_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("SAT", "DDSIM")})
        d_sat_clue = process_averages("q_sat_clue", remove_outliers=False, add_times=False, cpp=True)["red. ratio"].droplevel(0).to_frame().rename(columns={"red. ratio" : ("SAT", "d-CLUE")})
        for i in range(len(d_sat_clue)):
            d_sat_clue.iloc[i][0] = 2**d_sat_clue.iloc[i].name * d_sat_clue.iloc[i][0]
        d_sat_ddsim = process_averages("q_sat_ddsim", remove_outliers=False, add_times=False, cpp=True)["red. ratio"].droplevel(0).to_frame().rename(columns={"red. ratio" : ("SAT", "d-DD")})
        for i in range(len(d_sat_ddsim)):
            d_sat_ddsim.iloc[i][0] = 2**d_sat_ddsim.iloc[i].name * d_sat_ddsim.iloc[i][0]
        time_maxcut_clue = process_averages("q_maxcut_clue", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("MAXCUT", "CLUE")})
        time_maxcut_ddsim = process_averages("q_maxcut_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("MAXCUT", "DDSIM+CLUE")})
        time_maxcut_ddsim_alone = process_averages("q_maxcut_full_ddsim", remove_outliers=False, add_times=False, cpp=True)["time"].droplevel(0).to_frame().rename(columns={"time" : ("MAXCUT", "DDSIM")})
        d_maxcut_clue = process_averages("q_maxcut_clue", remove_outliers=False, add_times=False, cpp=True)["red. ratio"].droplevel(0).to_frame().rename(columns={"red. ratio" : ("MAXCUT", "d-CLUE")})
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

        data = reduce(
            lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
            columns)
    
    ## Now 'data' contains the table to be formatted
    data.index.name = None
    data.columns = pd.MultiIndex.from_tuples(data.columns, names=("",""))
    count_per_column = [len([c for c in data.columns if c[0] == l]) for l in data.columns.levels[0]]
    ## We generate the output
    if out == "text":
        if formt == "sci": # We use scientific notation
            pd.set_option("display.float_format", '{:.3E}'.format)
        elif formt == "float": # we print the float number
            pd.set_option("display.float_format", '{:.6f}'.format)
        print(data)
    else: # We need to format the output
        styler = data.style
        ## Place for formatting the table
        styler = styler.format({k: '{:.2f}' if k[1].startswith("d") or any(v == "d" for v in k) else '{:.3E}' for k in data.columns})
        
        ## Deciding the output of the table
        if out == "latex":
            with open(os.path.join(SCRIPT_DIR, "table_examples.tex"), "w") as file:
                styler.format_index("\\textbf{{{}}}", escape="latex", axis=1).to_latex(
                    file,
                    convert_css=True,
                    column_format="c|" + "|".join(e*"r" for e in count_per_column),
                    position="tp",
                    position_float="centering",
                    hrules=True,
                    label="tab:examples",
                    multirow_align="t",
                    multicol_align="c",
                    sparse_index=True,
                    sparse_columns=True,
                    caption=r"Evaluation of (single-step) quantum benchmarks from repository~\cite{quetschlich2022mqtbench}. The simulation \
                        times of DDSIM refer to the computation with respect to input $\ket{0}$, while the third and fourth columns report \
                        dimensions of bisimulation reductions. Instead, the fifth column reports the average computation time of $U \ket{x}$ \
                        via a reduction wrt $S_{\ket{x}}$, including the computation time of the bisimulation. A cumulative timeout of 500\,s \
                        is denoted by \textbf{TO}."
                )
        elif out == "html":
            with open(os.path.join(SCRIPT_DIR, "table_examples.html")) as file:
                styler.to_html(file)
        elif out == "style":
            return styler

def table2(*, cpp:bool=False, out:str="text", formt="sci"):
    ## Creating the table (depends on C++ argument)
    if not cpp:
        red_zero = process_averages("q_benchmark_clue", remove_outliers=False, observable=0)["red. ratio"].droplevel(-1).to_frame().rename(columns={"red. ratio": r"RR for $S_{\ket{0}}$"})
        avg_clue = process_averages("q_benchmark_clue", remove_outliers=False)[["red. ratio","time_lumping"]].rename(columns={"red. ratio": r"RR across $S_{\ket{x}}$", "time_lumping": "Avg. time (s)"})
        avg_time_ddsim = process_averages("q_benchmark_full_ddsim", remove_outliers=False, kappa=1)["time_iteration"].droplevel(-1).to_frame().rename(columns={"time_iteration": "DDSIM time"})

        data = reduce(
            lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
            [red_zero, avg_clue, avg_time_ddsim])    
        count_per_column = [4]
    else:
        red_zero_clue = process_averages("q_benchmark_clue", remove_outliers=False, observable=0, cpp=True)["red. ratio"].droplevel(-1).to_frame().rename(columns={"red. ratio": ("CLUE",r"RR for $S_{\ket{0}}$")})
        avg_clue_clue = process_averages("q_benchmark_clue", remove_outliers=False, cpp=True)[["red. ratio","time_lumping"]].rename(columns={"red. ratio": ("CLUE",r"RR across $S_{\ket{x}}$"), "time_lumping": ("CLUE", "Avg. time (s)")})
        # avg_clue_clue = process_averages("q_benchmark_clue", remove_outliers=False, cpp=True)[["time_lumping"]].rename(columns={"time_lumping": ("CLUE", "Avg. time (s)")})
        
        red_zero_ddsim = process_averages("q_benchmark_ddsim", remove_outliers=False, observable=0, cpp=True)["red. ratio"].droplevel(-1).to_frame().rename(columns={"red. ratio": ("DDSIM+CLUE",r"RR for $S_{\ket{0}}$")})
        avg_clue_ddsim = process_averages("q_benchmark_ddsim", remove_outliers=False, cpp=True)[["red. ratio","time_lumping"]].rename(columns={"red. ratio": ("DDSIM+CLUE",r"RR across $S_{\ket{x}}$"), "time_lumping": ("DDSIM+CLUE", "Avg. time (s)")})
        avg_time_ddsim = process_averages("q_benchmark_full_ddsim", remove_outliers=False, kappa=1, cpp=True)["time_iterations"].to_frame().rename(columns={"time_iterations": ("DDSIM", "Time")})

        data = reduce(
            lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
            [
                red_zero_clue, avg_clue_clue, 
                red_zero_ddsim, avg_clue_ddsim,
                avg_time_ddsim
            ])
        data.columns = pd.MultiIndex.from_tuples(data.columns, names=(None,None))
        count_per_column = [len([c for c in data.columns if c[0] == data.columns[0][0]])]
        while len(count_per_column) < len(data.columns.levels[0]):
            count_per_column += [len([c for c in data.columns if c[0] == data.columns[sum(count_per_column)][0]])]

    ## Now 'data' contains the table to be formatted
    data.index.names = [None, None]
    for column in data.columns:
        if "RR" in column or any("RR" in v for v in column):
            data[column] = data[column].apply(lambda x : 100*x)
    ## We generate the output
    if out == "text":
        if formt == "sci": # We use scientific notation
            pd.set_option("display.float_format", '{:.3E}'.format)
        elif formt == "float": # we print the float number
            pd.set_option("display.float_format", '{:.6f}'.format)
        print(data)
    else: # We need to format the output
        styler = data.style
        ## Place for formatting the table
        styler = styler.format({k: '{:.2f}\%' if "RR" in k or any("RR" in v for v in k) else '{:.3E}' for k in data.columns})
        if cpp: # we mark the minimum between the two avg. times
            styler = styler.highlight_min(subset=[("CLUE","Avg. time (s)"), ("DDSIM+CLUE","Avg. time (s)")], axis=1, props="font-weight:bold;")
        
        ## Deciding the output of the table
        if out == "latex":
            with open(os.path.join(SCRIPT_DIR, "table_benchmark.tex"), "w") as file:
                styler.format_index("\\textbf{{{}}}", axis=1).to_latex(
                    file,
                    convert_css=True,
                    column_format="rc|" + "|".join(e*"r" for e in count_per_column),
                    position="tp",
                    position_float="centering",
                    hrules=True,
                    clines="skip-last;data",
                    label="table:benchmark",
                    multirow_align="c",
                    multicol_align="c",
                    sparse_index=True,
                    sparse_columns=True,
                    caption=r"Comparison of simulation times between DDSIM and the reduced model by CLUE. The latter includes\
                              the runtimes for computing the bisimulations by Algorithm~\ref{alg}.  For SAT and MaxCut, column\
                              $d$ reports the average size of the reduced circuit and its theoretical bounds from Theorem~\ref{thm:ham:sat}, \
                              separated by backslash."
                )
        elif out == "html":
            with open(os.path.join(SCRIPT_DIR, "table_examples.html")) as file:
                styler.to_html(file)
        elif out == "style":
            return styler

def table3(*, cpp:bool=False, out:str="text", formt="sci"):
    ## Creating the table (depends on C++ argument)
    if not cpp:
        raise NotImplementedError
    else:
        zero_clue = process_averages("q_benchmark_clue", remove_outliers=False, observable='0', cpp=True).droplevel(-1)[["red. ratio", "time_lumping"]].rename(
            columns={"red. ratio": ("CLUE",r"RR for $S_{\ket{0}}$"), "time_lumping": ("CLUE", r"Time (s)")})
        zero_ddsim_clue = process_averages("q_benchmark_ddsim", remove_outliers=False, observable='0', cpp=True).droplevel(-1)[["red. ratio", "time_lumping"]].rename(
            columns={"red. ratio": ("DDSIM+CLUE",r"RR for $S_{\ket{0}}$"), "time_lumping": ("DDSIM+CLUE", r"Time (s)")})
        zero_ddsim_alone = process_averages("q_benchmark_full_ddsim", remove_outliers=False, observable=0, kappa=1, cpp=True).droplevel(-1)["time_iterations"].to_frame().rename(columns={"time_iterations": ("DDSIM", "Time (s)")})
        
        data = reduce(
            lambda p, q : p.merge(q, how="outer", left_index=True, right_index=True), 
            [
                zero_clue, zero_ddsim_clue, zero_ddsim_alone
            ])
        data.columns = pd.MultiIndex.from_tuples(data.columns, names=(None,None))
        count_per_column = [len([c for c in data.columns if c[0] == data.columns[0][0]])]
        while len(count_per_column) < len(data.columns.levels[0]):
            count_per_column += [len([c for c in data.columns if c[0] == data.columns[sum(count_per_column)][0]])]

    ## Now 'data' contains the table to be formatted
    data.index.names = [None, None]
    for column in data.columns:
        if "RR" in column or any("RR" in v for v in column):
            data[column] = data[column].apply(lambda x : 100*x)
    ## We generate the output
    if out == "text":
        if formt == "sci": # We use scientific notation
            pd.set_option("display.float_format", '{:.3E}'.format)
        elif formt == "float": # we print the float number
            pd.set_option("display.float_format", '{:.6f}'.format)
        print(data)
    else: # We need to format the output
        styler = data.style
        ## Place for formatting the table
        styler = styler.format({k: '{:.2f}\%' if "RR" in k or any("RR" in v for v in k) else '{:.3E}' for k in data.columns})
        if cpp: # we mark the minimum between the two avg. times
            styler = styler.highlight_min(subset=[("CLUE","Time (s)"), ("DDSIM+CLUE","Time (s)")], axis=1, props="font-weight:bold;")
        
        ## Deciding the output of the table
        if out == "latex":
            with open(os.path.join(SCRIPT_DIR, "table_benchmark_zeros.tex"), "w") as file:
                styler.format_index("\\textbf{{{}}}", axis=1).to_latex(
                    file,
                    convert_css=True,
                    column_format="rc|" + "|".join(e*"r" for e in count_per_column),
                    position="tp",
                    position_float="centering",
                    hrules=True,
                    clines="skip-last;data",
                    label="table:benchmark",
                    multirow_align="c",
                    multicol_align="c",
                    sparse_index=True,
                    sparse_columns=True,
                    caption=r"Some nice caption."
                )
        elif out == "html":
            with open(os.path.join(SCRIPT_DIR, "table_examples.html")) as file:
                styler.to_html(file)
        elif out == "style":
            return styler

#######################################################################################
### METHODS TO PROCESS THE ARGUMENTS OF THE SCRIPT
#######################################################################################
def avg_args(*argv):
    n = 0; observable = "all"; kappa = "all"; skip_kappa = None; remove_outliers = True; without_infinity = False; add_times = False; is_cpp = False
    while n < len(argv):
        if argv[n].startswith("-"):
            if argv[n].endswith("obs"):
                observable = argv[n+1]; n += 2
            elif argv[n].endswith("nk"):
                skip_kappa = int(argv[n+1]); n += 2
            elif argv[n].endswith("k"):
                kappa = int(argv[n+1]); n += 2
            elif argv[n].endswith("wo"):
                remove_outliers = False; n+=1
            elif argv[n].endswith("noinf"):
                without_infinity = True; n+=1
            elif argv[n].endswith("add"):
                add_times = True; n+=1
            elif argv[n].endswith("cpp"):
                is_cpp = True; n+=1
            else:
                n += 1
        else:
            n += 1
    return observable, kappa, skip_kappa, remove_outliers, without_infinity, add_times, is_cpp

def table_args(*argv):
    n=0; out="text"; formt="sci"; cpp=False
    while n < len(argv):
        if argv[n] == "-out":
            if argv[n+1] in ("text", "latex", "style", "html"):
                out = argv[n+1]
            n+=2
        elif argv[n] == "-format":
            if argv[n+1] in ("sci", "float"):
                formt = argv[n+1]
            n+=2
        elif argv[n] == "-cpp":
            cpp = True
            n+=1
        else:
            n+=1
    return out, formt, cpp

#######################################################################################
### MAIN SCRIPT
#######################################################################################  
if __name__ == "__main__":
    what = sys.argv[1]
    if what == "table1":
        out, formt, cpp = table_args(*sys.argv[2:])
        table1(out=out,formt=formt,cpp=cpp)
    elif what == "table2":
        out, formt, cpp = table_args(*sys.argv[2:])
        table2(out=out,formt=formt,cpp=cpp)
    elif what == "table3":
        out, formt, cpp = table_args(*sys.argv[2:])
        table3(out=out,formt=formt,cpp=cpp)
    else:
        observable, kappa, skip_kappa, remove_outliers, without_infinity, add_times, is_cpp = avg_args(*sys.argv[:2])

        print(
            process_averages(what, observable=observable, kappa=kappa, skip_kappa=skip_kappa, remove_outliers=remove_outliers, without_infinity=without_infinity, add_times=add_times,cpp=is_cpp)
        )
