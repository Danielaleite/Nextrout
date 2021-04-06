### import needed stuff



def forcing_generator(forcing_flag, grid, coord, extra_info):

	if forcing_flag == 'dirac':
		
		Nplus = extra_info['Nplus']
		Nminus = extra_info['Nminus']

		fplus = extra_info['fplus']
		fminus=  extra_info['fminus']

		xplus = extra_info['xplus'] # [[0.1,0.2],[0.3,0.4],[0.1,0.7]]
		xminus = extra_info['xminus'] # [[0.6,0.2],[0.8,0.4]]

		# set array forcing_dirac "evoluation" f=f^{+}-f^{-} on grid nodes
		forcing=np.zeros(grid.nnode)
		for i in range(Nplus):
		    inode=mt.Inode(coord,xplus[i])
		    forcing[inode]=fplus[i]
		for i in range(Nminus):
		    inode=mt.Inode(coord,xminus[i])
		    forcing[inode]=-fminus[i]

	elif forcing_flag == 'rect_cnst':

		rectangles_source = extra_info[0] # it should be a list of lists: every sublist has 3 elements: (x,y), w, h
		rectangles_sink = extra_info[1]
		
		forcing =np.zeros(grid.nnode)
		# iterate over the center of the triangles
		for cent in centers.keys():

		  x,y = centers[cent]

		  for rect in rectangles_source:

		    x_source = rect[0][0]
		    y_source = rect[0][1]
		    wo = rect[1]
		    ho = rect[2]

		    if (x_source <= x <= x_source+wo) and (y_source <= y <= y_source+ho):
		      inode=mt.Inode(coord,centers[cent])
		      #print(inode)
		      forcing[inode]=1

		  for rect in rectangles_sink:

		    x_sink = rect[0][0]
		    y_sink = rect[0][1]
		    wi = rect[1]
		    hi = rect[2]

		    if (x_sink <= x <= x_sink+wi) and (y_sink <= y <= y_sink+hi):
		      inode=mt.Inode(coord,centers[cent])
		      #print(inode)
		      forcing[inode]=-1
	else:
		raise ValueError('forcing not defined')

	return forcing

def kappa_generator(coord, kappa_flag):

	if kappa_flag == 'random':
	    x=coord[0]; y=coord[1]
	    fvalue=0.1+np.random.random_sample()
	else:
		raise ValueError('kappa not defined!')

    return fvalue


def dmk_cont(ndiv, forcing, beta, tdens0 = None,  nref= 0, flag_grid = 'unitsquare', kappa_flag = None):

	### building the grid
	# set mesh size 
	length=1.0/float(ndiv)

	# build grid using prebuild examples 
	points, vertices, coord,topol,element_attributes = example_grid.example_grid(flag_grid,length)
	# initialized fortran variable for the spatial discretization
	[grid,subgrid]=dmk_p1p0.init_geometry(topol, coord, 1)
	ncell=grid.ncell
	ntdens=grid.ncell
	npot=subgrid.nnode

	### defining dmk controls
	ctrl = Dmkcontrols.DmkCtrl()
	## redef param.
	ctrl.tolerance_system_variation = .000001
	print(ctrl.tolerance_system_variation)

	### defining input constructor
	dmkin=Dmkinputsdata.DmkInputs()
	Dmkinputsdata.dmkinputs_constructor(dmkin,0,100,100,True)
	## redef param.
	dmkin.pflux = beta
	print(dmkin.pflux)

	tdpot=Tdenspotentialsystem.tdpotsys()
	Tdenspotentialsystem.tdpotsys_constructor(tdpot,0,ntdens, npot,1)

	### defining tdens0

	if tdens0 is None:
		tdpot.tdens[:]=1.0
	else:
		tdpot.tdens = tdens0

	### defining kappa
	if kappa_flag is not None:
		# compute functions on cell centroids
		ncell=len(topol)
		bar_cell=mt.make_bar(coord,topol).transpose()
		kappa_cell=np.zeros([ncell]);
		for i in range(ncell):
		    kappa_cell[i] = kappa_generator(bar_cell[:,i])

	### defining rhs

	# initial integrated forcing term
	rhs=np.zeros(subgrid.ncell)
	dmk_p1p0.build_subgrid_rhs(subgrid,rhs, np.zeros(ncell),forcing)

	### init type for storing evolution/algorithm info
	timefun=Timefunctionals.evolfun()
	Timefunctionals.evolfun_constructor(timefun, 0,
	                                        ctrl.max_time_iterations,
	                                        ctrl.max_nonlinear_iterations)

	### solve with dmk
	info=0
	dmkp1p0_steady_data(grid, subgrid, tdpot, dmkin, ctrl, info,timefun=timefun)

	print(tdpot.tdens)

	return tdpot