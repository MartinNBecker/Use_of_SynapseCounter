//change dir to folder location with raw unstacked czi files
dir = "S:/shared/K6_Tammimies/CASK_project/Raw_data/IF_czi_pictures/Synapsin_Homer_pictures";

//change destination to folder where results should be saved
destination = "D:/Martin_KI/RunningProjects/CASK/ImmunoFluorescence/Puncta_quantification/Homer_Syn_Quantification/ResultFiles/";

//get list of files in image folder and only select czi files
list = newArray(0); 
files = getFileList(dir);

for(i=0;i<files.length;i++){ 
if(endsWith(files[i],".czi")){ 
        list=Array.concat(list,files[i]); 
        } 
} 


setBatchMode(true);

//for loop to process every image file one-by-one
for (k = 0; k < list.length; k++)
{
filename = list[k];
print(filename);
open(filename);
file=getTitle();
dotIndex = indexOf(file, " - "); 
file = substring(file, 0, dotIndex); 

//this part is analysing channel=0 and can be adapted to your prefered protocol
selectWindow(file + " - C=0"); //Red-Homer
run("Z Project...", "projection=[Max Intensity]"); //Max intensity of z-stacks
run("Enhance Contrast...", "saturated=0.35"); //enhance contrast to improve threshold
setAutoThreshold("Moments dark"); //sets threshold and ...
run("Convert to Mask"); //... creates new image with only 0 and 1 intensity values
run("Analyze Particles...", "size=1-Infinity show=Nothing display"); //count ALL particles. Inclusion criteria for particles will be set later in R
selectWindow("Results"); //select Results table window and ...
saveAs("Results", destination + filename + "_HomerResults.csv"); //... saves it with new extension
run("Close"); //closes the results table (!)
close(); //closes the image


//this part is analysing channel=1 and can be adapted to your prefered protocol
selectWindow(file + " - C=1"); //Green-SynapsinI
run("Z Project...", "projection=[Max Intensity]"); //Max intensity of z-stacks
run("Enhance Contrast...", "saturated=0.35"); //enhance contrast to improve threshold
setAutoThreshold("Moments dark"); //sets threshold and ...
run("Convert to Mask"); //... creates new image with only 0 and 1 intensity values
run("Analyze Particles...", "size=1-Infinity show=Nothing display"); //count ALL particles. Inclusion criteria for particles will be set later in R
selectWindow("Results"); //select Results table window and ...
saveAs("Results", destination + filename + "_SynapsinResults.csv"); //... saves it with new extension
run("Close"); //closes the results table (!)
close(); //closes the image


//this part is analysing channel=2 and can be adapted to your prefered protocol
selectWindow(file + " - C=2"); //DAPI
run("Z Project...", "projection=[Max Intensity]");
run("Enhance Contrast...", "saturated=0.35");
setAutoThreshold("Triangle dark");
run("Convert to Mask");
run("Analyze Particles...", "size=1-Infinity show=Nothing display");
selectWindow("Results"); 
saveAs("Results", destination + filename + "_NucleiResults.csv");
run("Close");
close();
}

print("Finished");
setBatchMode(false);
