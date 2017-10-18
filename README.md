!To run SCIP [java version 8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) is needed!

If you have an older version follow [this tutorial](https://de.mathworks.com/matlabcentral/answers/130359-how-do-i-change-the-java-virtual-machine-jvm-that-matlab-is-using-on-windows).

1. Replace the inputImage variable (first line) with the Image to be processed. 
2. Run SCIP.m

All processed steps will be displayed in the command window.

Outputs will be stored in the Results_'Imagename' folder

Outputfiles:

- Coordinates.xslx stores all x,y,z coordinates of detected cells.
- MAX_proj.tif will show you all detected cells on the projection.
- stackCropped.tif is the cropped original stack using the envelope around the fitted surface