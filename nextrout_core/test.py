import os
from main import nextrout

### testing:

# generate forcing

forcing_flag = "dirac"

if forcing_flag == "rect_cnst":
    x_source1, y_source1 = (0.2, 0.2)
    x_source2, y_source2 = (0.2, 0.7)
    wo1 = 0.05
    ho1 = 0.1
    rectangles_source = [
        [(x_source1, y_source1), wo1, ho1],
        [(x_source2, y_source2), wo1, ho1],
    ]  # bottom left cornner, width, height

    x_sink, y_sink = (0.8, 0.8)
    wi = 0.1
    hi = 0.1
    rectangles_sink = [[(x_sink, y_sink), wi, hi]]

    extra_info = [rectangles_source, rectangles_sink]

elif forcing_flag == "rect_cnst2":
    x_source1, y_source1 = (0.1, 0.1)
    wo1 = 0.2
    ho1 = 0.5
    rectangles_source = [
        [(x_source1, y_source1), wo1, ho1]
    ]  # bottom left cornner, width, height

    x_sink, y_sink = (0.7, 0.1)
    wi = 0.2
    hi = 0.5
    rectangles_sink = [[(x_sink, y_sink), wi, hi]]

    extra_info = [rectangles_source, rectangles_sink]

elif forcing_flag == "rect_cnst_d":
    x_source1, y_source1 = (0.1, 0.1)
    wo1 = 0.05
    ho1 = 0.05
    rectangles_source = [
        [(x_source1, y_source1), wo1, ho1]
    ]  # bottom left cornner, width, height

    x_sink, y_sink = (0.7, 0.1)
    wi = 0.05
    hi = 0.05
    rectangles_sink = [[(x_sink, y_sink), wi, hi]]

    extra_info = [rectangles_source, rectangles_sink]


elif forcing_flag == "dirac":
    Nplus = 3
    Nminus = 2

    fplus = [1, 2, 3]
    fminus = [4, 2]

    xplus = [[0.1, 0.21], [0.3, 0.4], [0.1, 0.7]]
    xminus = [[0.6, 0.2], [0.8, 0.4]]

    extra_info = {
        "Nplus": Nplus,
        "Nminus": Nminus,
        "fplus": fplus,
        "fminus": fminus,
        "xplus": xplus,
        "xminus": xminus,
    }

elif forcing_flag == "dirac2":
    fplus = [1, 2, 3, 1]
    fminus = [4, 2, 1]

    xplus = [[0.1, 0.2], [0.3, 0.4], [0.1, 0.7], [0.1, 0.9]]
    xminus = [[0.6, 0.2], [0.8, 0.4], [0.9, 0.5]]

    Nplus = len(xplus)
    Nminus = len(xminus)

    extra_info = {
        "Nplus": Nplus,
        "Nminus": Nminus,
        "fplus": fplus,
        "fminus": fminus,
        "xplus": xplus,
        "xminus": xminus,
    }

elif forcing_flag == "dirac3":
    fplus = [1]
    fminus = [1]

    xplus = [[0.1, 0.2]]
    xminus = [[0.6, 0.2]]

    Nplus = len(xplus)
    Nminus = len(xminus)

    extra_info = {
        "Nplus": Nplus,
        "Nminus": Nminus,
        "fplus": fplus,
        "fminus": fminus,
        "xplus": xplus,
        "xminus": xminus,
    }

elif forcing_flag == "dirac4":
    fplus = [3, 2]
    fminus = [4, 2]

    xplus = [[0.1, 0.21], [0.3, 0.9]]
    xminus = [[0.6, 0.2], [0.8, 0.4]]

    Nplus = len(xplus)
    Nminus = len(xminus)

    extra_info = {
        "Nplus": Nplus,
        "Nminus": Nminus,
        "fplus": fplus,
        "fminus": fminus,
        "xplus": xplus,
        "xminus": xminus,
    }


flags = [
    "whole_convex_hull+btns_centr",
    "branch_convex_hull+btns_centr",
    "btns_centr",
    "single",
]
test_flag = forcing_flag
### running nextrout

if test_flag == "dirac":
    beta_c = 1.5
    beta_d = 1.5
    ndiv = 20
    graph_type = "1"
    weighting_method = "ER"
    min_pe = 1e-2
    min_f = 1e-2
    BPw = "tdens"
    weighting_method_simplification = "BPW"
    stop_thresh_f = 1e-12
    verbose = False
    weight_flag = "length"
    btns_factor_source = 1
    btns_factor_sink = 1
    terminal_criterion = flags[3]
    storing = "outputs/"

elif test_flag == "dirac2":
    beta_c = 1.2
    beta_d = 1.8
    ndiv = 20
    graph_type = "1"
    weighting_method = "ER"
    min_pe = 1e-2
    min_f = 1e-2
    BPw = "tdens"
    weighting_method_simplification = "BPW"
    stop_thresh_f = 1e-12
    verbose = False
    weight_flag = "length"
    btns_factor_source = 1
    btns_factor_sink = 1
    terminal_criterion = flags[3]
    storing = "outputs/"

elif test_flag == "dirac3":
    beta_c = 1.5
    beta_d = 1.5
    ndiv = 15
    graph_type = "1"
    weighting_method = "ER"
    min_pe = 1e-2
    min_f = 1e-2
    BPw = "flux"
    weighting_method_simplification = "BPW"
    stop_thresh_f = 1e-12
    verbose = False
    weight_flag = "length"
    btns_factor_source = 1
    btns_factor_sink = 1
    terminal_criterion = flags[3]
    storing = "outputs/"


os.system("mkdir " + storing + test_flag)
storing = storing + "/" + test_flag

nextrout(
    forcing_flag,
    extra_info,
    beta_c=beta_c,
    beta_d=beta_d,
    ndiv=ndiv,
    graph_type=graph_type,
    weighting_method=weighting_method,
    min_pe=min_pe,
    min_f=min_f,
    BPw=BPw,
    weighting_method_simplification=weighting_method_simplification,
    stop_thresh_f=stop_thresh_f,
    verbose=verbose,
    weight_flag=weight_flag,
    btns_factor_source=btns_factor_source,
    btns_factor_sink=btns_factor_sink,
    terminal_criterion=terminal_criterion,
    storing=storing,
)

print("Nextrout algorithm finished!")
