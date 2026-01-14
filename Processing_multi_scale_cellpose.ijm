/*
 * Macro template to process multiple images in a folder
 */

// input parameters
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ Integer (label = "Droplet channnel", value = 1, style = "spinner") droplet_channel
#@ File (label = "Input Cellpose environment path", style = "directory") env_path
#@ Integer (label = "Minimum size for droplets (pixels)", value = 10, style = "spinner") min_size_droplet
#@ Integer (label = "Cell size for Cellpose at small scale", value = 5, style = "spinner") cell_scale_small
#@ Float (label = "Threshold for Cellpose at small scale", value = 0.0, style="format:#.##") cell_threshold_small
#@ Integer (label = "Cell size for Cellpose at medium scale", value = 20, style = "spinner") cell_scale_medium
#@ Float (label = "Threshold for Cellpose at medium scale", value = 0.0, style="format:#.##") cell_threshold_medium
#@ Integer (label = "Cell size for Cellpose at large scale", value = 50, style = "spinner") cell_scale_large
#@ Float (label = "Threshold for Cellpose at large scale", value = 0.0, style="format:#.##") cell_threshold_large
#@ String (label = "File suffix", value = ".tif") suffix

// call to the main function "processFolder"
processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	///////////// initial cleaning /////////////////
	// close all images
	run("Close All");
	// reset ROI manager
	roiManager("Reset");
	// clear results
	run("Clear Results");
	if ( cell_scale_small > 0 ) {

		///////////// apply pipeline to input images /////////////////
		// get the files in the input folder
		list = getFileList(input);
		list = Array.sort(list);
		// loop over the files
		for (i = 0; i < list.length; i++) {
			// if current file ends with the suffix given as input parameter, call function "processFile" to process it
			if(endsWith(list[i], suffix))
				processFile(input, output, list[i]);
		}
		// reset ROI manager
		roiManager("Reset");
		// save parameters
		// create results table
		Table.create("Results");
		setResult("Droplet channnel", 0, droplet_channel);
		setResult("Minimum size for droplets (pixels)", 0, min_size_droplet);
		setResult("Cell size for Cellpose at small scale", 0, cell_scale_small);
		setResult("Threshold for Cellpose at small scale", 0, cell_threshold_small);
		setResult("Cell size for Cellpose at medium scale", 0, cell_scale_medium);
		setResult("Threshold for Cellpose at medium scale", 0, cell_threshold_medium);
		setResult("Cell size for Cellpose at large scale", 0, cell_scale_large);
		setResult("Threshold for Cellpose at large scale", 0, cell_threshold_large);
		updateResults();
		// save results
		saveAs("Results", output + File.separator + "parameters.csv");	
		// close results table
		selectWindow("Results"); 
		run("Close");
	}
	else{
		IJ.log("Cell size for Cellpose at small scale must be superior to 0");
	}

}

function processFile(input, output, file) {

	// open image
	open(input + File.separator + file);
	rename("wholeInput");
	// split channels
	run("Split Channels");
	// select droplet channel
	selectImage("C" + droplet_channel + "-wholeInput");
	// rename
	rename("InputVolume");
	// get dimensions
	getDimensions(width, height, channels, nbSlices, frames);
	// get resolution
	getPixelSize(scaleUnit, pixelWidth, pixelHeight);
	
	// create results table
	Table.create("Results");
	for (z = 1; z <= nbSlices; z++) {
		// dulicate current slice
		selectImage("InputVolume");
		setSlice(z);
		run("Duplicate...", "title=Input");	
		// increase contrast
		run("Duplicate...", "title=InputForVisualInspection");	
		run("Enhance Contrast", "saturated=0.35");
		run("Apply LUT");
		
		current_nb_droplets = 0;
		// create new image with regions to exclude
		newImage("RegionsToExclude", "32-bit black", width, height, 1);
		run("Add...", "value=1");
	
		// run cellpose at large scale
		selectImage("Input");
		if (cell_scale_large > (-1) ) {
			// segment cells with cellpose
			run("Cellpose ...", "env_path=" + env_path + " env_type=venv model=cyto2 model_path=path\\to\\own_cellpose_model diameter=" + cell_scale_large + " ch1=0 ch2=-1 additional_flags=[--use_gpu, --cellprob_threshold, " + cell_threshold_large + ", --flow_threshold, 2.0]");
			// test if any droplet was identified
			if ( isOpen("Input-cellpose") ) {
				selectImage("Input-cellpose");
				rename("LargeDroplets");
				// get droplets at current scale as a binary mask
				setThreshold(1, 65535, "raw");
				run("Create Mask");
				// invert it to add to regions to exclude at other scales
				run("Invert");
				run("Divide...", "value=255");
				rename("LargeDropletsToExclude");
				imageCalculator("Multiply create 32-bit", "LargeDropletsToExclude","RegionsToExclude");
				close("RegionsToExclude");	
				rename("RegionsToExclude");	
				// close useless images
				close("Input-cellpose"); 
				close("LargeDropletsToExclude"); 
				//selectImage("LargeDroplets");
				//getStatistics(first, intensity, min, current_nb_droplets, std, histogram);
			}
		}
		// run LoG at medium scale
		selectImage("Input");
		if (cell_scale_medium > (-1) ) {
			// segment cells with cellpose
			run("Cellpose ...", "env_path=" + env_path + " env_type=venv model=cyto2 model_path=path\\to\\own_cellpose_model diameter=" + cell_scale_medium + " ch1=0 ch2=-1 additional_flags=[--use_gpu, --cellprob_threshold, " + cell_threshold_medium + ", --flow_threshold, 2.0]");
			// test if any droplet was identified
			if ( isOpen("Input-cellpose") ) {
				// remove saturated areas
				imageCalculator("Multiply create 32-bit", "Input-cellpose","RegionsToExclude");
				// connected components
				run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
				rename("MediumDropletsBeforeIdUpdate1");
				run("Duplicate...", "title=MediumDropletsBeforeIdUpdate2");
				run("Add...", "value=" + current_nb_droplets + "");
				selectImage("MediumDropletsBeforeIdUpdate1");
				setThreshold(1, 65535, "raw");
				run("Create Mask");
				run("Divide...", "value=255");
				rename("MediumDropletsBeforeIdUpdate3");
				imageCalculator("Multiply create 32-bit", "MediumDropletsBeforeIdUpdate2","MediumDropletsBeforeIdUpdate3");
				rename("MediumDroplets");
				// get droplets at current scale as a binary mask
				setThreshold(1, 65535, "raw");
				run("Create Mask");
				// invert it to add to regions to exclude at other scales
				run("Invert");
				run("Divide...", "value=255");
				rename("MediumDropletsToExclude");
				imageCalculator("Multiply create 32-bit", "MediumDropletsToExclude","RegionsToExclude");
				close("RegionsToExclude");	
				rename("RegionsToExclude");	
				// close useless images
				close("Input-cellpose"); 
				close("Result of Input-cellpose"); 
				close("MediumDropletsToExclude");
				close("MediumDropletsBeforeIdUpdate1");
				close("MediumDropletsBeforeIdUpdate2");
				close("MediumDropletsBeforeIdUpdate3");
				//selectImage("MediumDroplets");
				//getStatistics(first, intensity, min, current_nb_droplets, std, histogram);
			}
		}
		// run LoG at small scale
		selectImage("Input");
		if (cell_scale_small > (-1) ) {
			// segment cells with cellpose
			run("Cellpose ...", "env_path=" + env_path + " env_type=venv model=cyto2 model_path=path\\to\\own_cellpose_model diameter=" + cell_scale_small + " ch1=0 ch2=-1 additional_flags=[--use_gpu, --cellprob_threshold, " + cell_threshold_small + ", --flow_threshold, 2.0]");
			// test if any droplet was identified
			if ( isOpen("Input-cellpose") ) {
				// remove saturated areas
				imageCalculator("Multiply create 32-bit", "Input-cellpose","RegionsToExclude");
				// connected components
				run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
				run("Label Size Filtering", "operation=Greater_Than size=" + min_size_droplet + "");
				rename("SmallDropletsBeforeIdUpdate1");
				run("Duplicate...", "title=SmallDropletsBeforeIdUpdate2");
				run("Add...", "value=" + current_nb_droplets + "");
				selectImage("SmallDropletsBeforeIdUpdate1");
				setThreshold(1, 65535, "raw");
				run("Create Mask");
				run("Divide...", "value=255");
				rename("SmallDropletsBeforeIdUpdate3");
				imageCalculator("Multiply create 32-bit", "SmallDropletsBeforeIdUpdate2","SmallDropletsBeforeIdUpdate3");
				rename("SmallDroplets");
				// get droplets at current scale as a binary mask
				setThreshold(1, 65535, "raw");
				run("Create Mask");
				// invert it to add to regions to exclude at other scales
				run("Invert");
				run("Divide...", "value=255");
				rename("SmallDropletsToExclude");
				imageCalculator("Multiply create 32-bit", "SmallDropletsToExclude","RegionsToExclude");
				close("RegionsToExclude");	
				rename("RegionsToExclude");	
				// close useless images
				close("Input-cellpose"); 
				close("Result of Input-cellpose"); 
				close("SmallDropletsToExclude");
				close("SmallDropletsBeforeIdUpdate1");
				close("SmallDropletsBeforeIdUpdate2");
				close("SmallDropletsBeforeIdUpdate3");
				//selectImage("SmallDroplets");
				//getStatistics(first, intensity, min, current_nb_droplets, std, histogram);
			}
		}
		
		setResult("Slice", z-1, z);
		updateResults();
		if ( isOpen("LargeDroplets") ) {
			// compute number of large droplets
			selectImage("LargeDroplets");
			run("Label image to ROIs", "rm=[RoiManager[]]");
			// add number to results
			setResult("#large droplets", z-1, RoiManager.size);
			updateResults();
			// measure average intensity over large droplets
			// threshold droplets
			selectImage("LargeDroplets");
			setThreshold(1.0000, 1000000000000000000000000000000.0000);
			// create ROI and add it to ROI manager
			run("Create Selection");
			roiManager("Add");
			// select corresponding input
			selectImage("Input");
			roiManager("Select", RoiManager.size-1);
			// measure average intensity
			getStatistics(large_droplets_area, large_droplets_avg_intensity, min, max, std, histogram);
			// add avg intensity and total area to results
			setResult("Average intensity for large droplets", z-1, large_droplets_avg_intensity);
			setResult("Total area for large droplets", z-1, large_droplets_area);
			updateResults();
			// remove ROI corresponding to all droplets
			roiManager("Select", RoiManager.size-1);
			roiManager("Delete");
			// change ROI color
			count = roiManager("count");
			array = newArray(count);
  			for (i=0; i<array.length; i++) {
      			array[i] = i;
  			}
			roiManager("select", array);
			roiManager("Set Color", "cyan");
			// add ROIs for visual inspection
			selectImage("InputForVisualInspection");
			roiManager("Show All without labels");
			run("Flatten");
			close("InputForVisualInspection");
			rename("InputForVisualInspection");
			close("LargeDroplets");
			// reset ROI manager
			roiManager("Reset");
		}
		if ( isOpen("MediumDroplets") ) {
			// compute number of medium droplets
			selectImage("MediumDroplets");
			run("Label image to ROIs", "rm=[RoiManager[]]");
			// add number to results
			setResult("#medium droplets", z-1, RoiManager.size);
			updateResults();
			// measure average intensity over medium droplets
			// threshold droplets
			selectImage("MediumDroplets");
			setThreshold(1.0000, 1000000000000000000000000000000.0000);
			// create ROI and add it to ROI manager
			run("Create Selection");
			roiManager("Add");
			// select corresponding input
			selectImage("Input");
			roiManager("Select", RoiManager.size-1);
			// measure average intensity
			getStatistics(medium_droplets_area, medium_droplets_avg_intensity, min, max, std, histogram);
			// add avg intensity and total area to results
			setResult("Average intensity for medium droplets", z-1, medium_droplets_avg_intensity);
			setResult("Total area for medium droplets", z-1, medium_droplets_area);
			updateResults();
			// remove ROI corresponding to all droplets
			roiManager("Select", RoiManager.size-1);
			roiManager("Delete");
			// change ROI color
			count = roiManager("count");
			array = newArray(count);
  			for (i=0; i<array.length; i++) {
      			array[i] = i;
  			}
			roiManager("select", array);
			roiManager("Set Color", "yellow");
			// add ROIs for visual inspection
			selectImage("InputForVisualInspection");
			roiManager("Show All without labels");
			run("Flatten");
			close("InputForVisualInspection");
			rename("InputForVisualInspection");
			close("MediumDroplets");
			// reset ROI manager
			roiManager("Reset");
		}
		if ( isOpen("SmallDroplets") ) {
			// compute number of small droplets
			selectImage("SmallDroplets");
			run("Label image to ROIs", "rm=[RoiManager[]]");
			// add number to results
			setResult("#small droplets", z-1, RoiManager.size);
			updateResults();
			// measure average intensity over small droplets
			// threshold droplets
			selectImage("SmallDroplets");
			setThreshold(1.0000, 1000000000000000000000000000000.0000);
			// create ROI and add it to ROI manager
			run("Create Selection");
			roiManager("Add");
			// select corresponding input
			selectImage("Input");
			roiManager("Select", RoiManager.size-1);
			// measure average intensity
			getStatistics(small_droplets_area, small_droplets_avg_intensity, min, max, std, histogram);
			// add avg intensity and total area to results
			setResult("Average intensity for small droplets", z-1, small_droplets_avg_intensity);
			setResult("Total area for small droplets", z-1, small_droplets_area);
			updateResults();
			// remove ROI corresponding to all droplets
			roiManager("Select", RoiManager.size-1);
			roiManager("Delete");
			// change ROI color
			count = roiManager("count");
			array = newArray(count);
  			for (i=0; i<array.length; i++) {
      			array[i] = i;
  			}
			roiManager("select", array);
			roiManager("Set Color", "red");
			// add ROIs for visual inspection
			selectImage("InputForVisualInspection");
			roiManager("Show All without labels");
			run("Flatten");
			close("InputForVisualInspection");
			rename("InputForVisualInspection");
			close("SmallDroplets");
			// reset ROI manager
			roiManager("Reset");
		}
		if ( isOpen("VisualInspectionOutput") ) {
			run("Concatenate...", "open image1=VisualInspectionOutput image2=InputForVisualInspection image3=[-- None --]");
			rename("VisualInspectionOutput");
			close("InputForVisualInspection");
		}
		else {
			selectImage("InputForVisualInspection");
			rename("VisualInspectionOutput");
		}
		close("Input");
		close("RegionsToExclude");
	}


	// save visual inspection
	if ( isOpen("VisualInspectionOutput") ) {
		selectImage("VisualInspectionOutput");
		file_name1 = replace(file, suffix, ".tif");
		saveAs("tif", output + File.separator + file_name1);	
	}
	
	// save results
	file_name2 = replace(file, suffix, "_results.csv");
	saveAs("Results", output + File.separator + file_name2);	

		
	///////////// clear everything /////////////////
	// close all images
	run("Close All");
	// reset ROI manager
	roiManager("Reset");
	// close results table
	close("Results"); 
}
