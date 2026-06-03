TestMethod        = True
SampleWindow      = 0      ## whole trajectory - single molecules
FirstDimension    = 0      ## 1-2 torsion
SecondDimension   = 1      ## 2-3 torsion
num_processors    = 1      ## 1 core
minima_depth      = 10     ## 10 minimas
windows           = {8000:'100', 7000:'150', 6000:'200', 5000:'250',0:'500'}
SampleWindowTime  = windows.get(SampleWindow)