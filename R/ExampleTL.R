#Loading tcx data
library(trackeR)
run2 = readTCX('ExampleData/activity_8377782029.tcx');

#Loading cwa axivity data
library(GGIRread)
ax3 = readAxivity('ExampleData/ax3.cwa', start=0, end=1000000)


