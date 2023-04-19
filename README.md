# trainingLoad
Installing the package requires the devtools package. Use the following command to install the package:

install.packages("devtools")


<B>Installing the trainingLoad package</B>
<br>
The R package can be installed using the following command

devtools::install_github("jbrond/trainingLoad")

<B>Example data usage</B>
<br>
Use the trackeR package to read TCX files:<br>

library(trackeR)<br>
garmin = readTCX('ExampleData/activity_8377782029.tcx')<br>

For loading Axivity Accelerometry files use the GGIRread package

library(GGIRread)<br>
ax3 = readAxivity('ExampleData/ax3.cwa', start=0, end=200)<br>

Start and end with readAxivity is in blocks of 3 seconds. Thus with end=200 and sampling frequency of 100 Hz we get the first 10 minues of data (200&#42;3&#42;100 is 60000 samples)
