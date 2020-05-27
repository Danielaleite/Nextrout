# Visit 2.10.0 log file
ScriptVersion = "2.10.0"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
RestoreSession("../../../runs/sessions/pflux300_overlap.session", 0)
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "/home/enrico/bbmuffe/p1_muffe/runs/figures"
SaveWindowAtts.fileName = "optdens_pflux300_overlap"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1600
SaveWindowAtts.height = 1200
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.advancedMultiWindowSave = 0
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
