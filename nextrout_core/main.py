#### import the needed stuff









def nextrout():

	# run the dmk_cont

	forcing = forcing_generator(forcing_flag, grid, coord, extra_info)
	tdpot = dmk_cont(ndiv, forcing, beta, tdens0,  nref= 0, flag_grid = 'unitsquare')

	# run the graph extraction

	Gpe = pre_extraction(coord, topol, tdpot, min_)

	# run the dmk_discrete

	Gf = filtering()

	# post process the output


	return Gf


### testing:

x_source1, y_source1 = (.2,.2)
x_source2, y_source2 = (.2,.7)
wo1 = .05
ho1 = .1
rectangles_source = [[(x_source1,y_source1),wo1,ho1],[(x_source2,y_source2),wo1,ho1]] # bottom left cornner, width, height

x_sink, y_sink = (.8,.8)
wi = .1
hi = .1
rectangles_sink = [[(x_sink,y_sink),wi,hi]] # bottom left cornner, width, height
