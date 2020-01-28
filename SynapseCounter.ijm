//Make sure to enable Windowless Importer for image files from BioFormats Plugin
//Install SynapseCounter plugin: https://github.com/SynPuCo/SynapseCounter
//Here SynapseCounter does not work headless and it is required to click the "OK" button for every image


//change dir to folder location with raw unstacked czi files
dir = "/path/to/images";

//change destination to folder where results should be saved
destination = "/path/to/ResultFiles/";

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

run("Z Project...", "projection=[Max Intensity]");
run("Synapse Counter");
selectWindow("SynapseCounter results"); 
saveAs("Results", destination + filename + "_SynapseCounterResults.csv");
run("Close");
close();

}

print("Finished");
setBatchMode(false);
