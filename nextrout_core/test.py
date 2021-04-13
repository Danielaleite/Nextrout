from main import nextrout
### testing:

# generate forcing

forcing_flag = 'rect_cnst'

if forcing_flag == 'rect_cnst':

    x_source1, y_source1 = (.2,.2)
    x_source2, y_source2 = (.2,.7)
    wo1 = .05
    ho1 = .1
    rectangles_source = [[(x_source1,y_source1),wo1,ho1],[(x_source2,y_source2),wo1,ho1]] # bottom left cornner, width, height

    x_sink, y_sink = (.8,.8)
    wi = .1
    hi = .1
    rectangles_sink = [[(x_sink,y_sink),wi,hi]]

    extra_info = [rectangles_source,rectangles_sink]

elif forcing_flag == 'rect_cnst2':

    x_source1, y_source1 = (.1,.1)
    wo1 = .2
    ho1 = .5
    rectangles_source = [[(x_source1,y_source1),wo1,ho1]] # bottom left cornner, width, height

    x_sink, y_sink = (.7,.1)
    wi = .2
    hi = .5
    rectangles_sink = [[(x_sink,y_sink),wi,hi]]

    extra_info = [rectangles_source,rectangles_sink]

elif forcing_flag == 'rect_cnst_d':

    x_source1, y_source1 = (.1,.1)
    wo1 = .05
    ho1 = .05
    rectangles_source = [[(x_source1,y_source1),wo1,ho1]] # bottom left cornner, width, height

    x_sink, y_sink = (.7,.1)
    wi = .05
    hi = .05    
    rectangles_sink = [[(x_sink,y_sink),wi,hi]]

    extra_info = [rectangles_source,rectangles_sink]

    

elif forcing_flag == 'dirac':

    Nplus=3
    Nminus=2

    fplus=[1,2,3]
    fminus=[4,2]

    xplus=[[0.1,0.21],[0.3,0.4],[0.1,0.7]]
    xminus=[[0.6,0.2],[0.8,0.4]]

    extra_info = {'Nplus':Nplus,
                   'Nminus':Nminus,
                    'fplus':fplus,
                    'fminus':fminus,
                    'xplus':xplus,
                    'xminus':xminus}

elif forcing_flag == 'dirac2':

    fplus=[1,2,3,1]
    fminus=[4,2,1]

    xplus=[[0.1,0.2],[0.3,0.4],[0.1,0.7], [0.1,0.9]]
    xminus=[[0.6,0.2],[0.8,0.4],[0.9,0.5]]

    Nplus = len(xplus)
    Nminus = len(xminus)

    extra_info = {'Nplus':Nplus,
                   'Nminus':Nminus,
                    'fplus':fplus,
                    'fminus':fminus,
                    'xplus':xplus,
                    'xminus':xminus}

elif forcing_flag == 'dirac3':

    fplus=[3,2]
    fminus=[4,2]

    xplus=[[0.1,0.21],[0.3,0.9]]
    xminus=[[0.6,0.2],[0.8,0.4]]

    Nplus = len(xplus)
    Nminus = len(xminus)

    extra_info = {'Nplus':Nplus,
                   'Nminus':Nminus,
                    'fplus':fplus,
                    'fminus':fminus,
                    'xplus':xplus,
                    'xminus':xminus}

beta_c = 1.5
beta_d = 1.5
flags = ['whole_convex_hull+btns_centr','branch_convex_hull+btns_centr','btns_centr','single']

### running nextrout

nextrout(forcing_flag,
    extra_info,
    beta_c,
    beta_d = beta_d, 
    ndiv = 30, 
    graph_type='1',
    weighting_method = 'ER',
    min_pe = 0.01,
    min_f = 0.1,
    BPw = 'flux',
    weighting_method_simplification ='ER',
    stop_thresh_f = 1e-8,
    verbose = False,
    weight_flag = 'length',
    btns_factor_source=.5, 
    btns_factor_sink=.5,
    terminal_criterion =  flags[3],
    storing = './outputs/')
