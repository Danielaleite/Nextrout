import sys

session_file=sys.argv[1]
print session_file

RestoreSession(str(session_file),0)
 
# Get the operator attributes for each operator and print them
for opIndex in range(len(GetPlotList().GetPlots(0).operators)):
    atts = GetOperatorOptions(opIndex)
    print atts

SaveWindow()

sys.exit()

