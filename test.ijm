// select cell 23 jan 26
run("Duplicate...", " ");
setBackgroundColor(10, 10, 10);//to the background camera
setForegroundColor(255, 255, 255);
run("Clear Outside");
//run("Select None");


run("Tubeness", "sigma=2");// change size of filter here
//run("Threshold...");
setAutoThreshold("Otsu dark");
run("Create Selection");
close();
run("Restore Selection");
run("Set Measurements...", "area mean integrated display redirect=None decimal=3");
run("Measure");

//run("Channels Tool...");
Property.set("CompositeProjection", "null");
Stack.setDisplayMode("grayscale");
Stack.setChannel(4);
//setTool("freehand");
run("Auto Threshold", "method=Otsu white");
run("Set Measurements...", "area mean standard integrated median display redirect=None decimal=3");
run("Measure");


//run("Channels Tool...");
Property.set("CompositeProjection", "null");
Stack.setDisplayMode("grayscale");
Stack.setChannel(4);


setOption("BlackBackground", false);
run("Convert to Mask", "method=Otsu background=Dark calculate");