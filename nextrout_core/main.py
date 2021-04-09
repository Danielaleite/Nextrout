#### import the needed stuff
import dmk_cont
import pre_extraction
import filtering

import matplotlib.pyplot as plt
import networkx as nx






def nextrout(forcing_flag,extra_info,beta,ndiv=20, verbose = True):

	# run the dmk_cont
	grid, subgrid, points, vertices, coord,topol,element_attributes = dmk_cont.grid_gen(ndiv)
	forcing = dmk_cont.forcing_generator(forcing_flag, grid, coord, topol, extra_info=extra_info)
	tdpot = dmk_cont.dmk_cont(ndiv,forcing,beta)
	print(tdpot.tdens)
	# run the graph extraction
	print(len(coord),len(topol))
	Gpe = pre_extraction.pre_extr(coord, topol, tdpot, min_= 0.00001)

	if verbose:
		pos = nx.get_node_attributes(Gpe,'pos')
		nx.draw(Gpe,pos, node_size = 10, node_color = 'k')
		plt.show()
	# run the dmk_discrete

	Gf = filtering.filtering(Gpe)
	if verbose:
		pos = nx.get_node_attributes(Gf,'pos')
		nx.draw(Gf,pos, node_size = 10, node_color = 'k')
		plt.show()

	# post process the output


	return Gf


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

beta = 1.5
nextrout(forcing_flag,extra_info,beta,ndiv=20 )
