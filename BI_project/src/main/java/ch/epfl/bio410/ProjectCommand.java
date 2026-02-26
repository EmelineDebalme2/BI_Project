package ch.epfl.bio410;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.Plot;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.ZProjector;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import net.imagej.ImageJ;
//import org.jfree.chart.title.LegendTitle;
import org.scijava.command.Command;
import org.scijava.plugin.Plugin;
import ij.gui.GenericDialog;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * This plugin aims at quantifying over time the Mycobacterium tuberculosis infection in macrophages.
 * The data is composed of 3 images: wild type 1, wild type 2, and mutant cells.
 * These images have slices with 2 µm z-step, frames with 2h frame interval, and 2 channels (one for macrophages,
 * one for bacteria).
 * It implements the following pipeline:
 * - Open the images through a dialog box
 * - Process and segment the images
 * - Do 4 types of measurements on the images to quantify the infection
 * - Plot the measurements
 */
@Plugin(type = Command.class, menuPath = "Plugins>BII>Infection Quantification")
public class ProjectCommand implements Command {

	public void run() {

		// Get the data path via a dialog box
		GenericDialog dlg = new GenericDialog("Get the data path");
		dlg.addDirectoryField("Data path", Prefs.getHomeDir());
		dlg.showDialog();
		if (dlg.wasCanceled()) return;
		String path = dlg.getNextString();

		//1. Open the images
		// Wild Type 1
		ImagePlus wt_1 = IJ.openImage(path + "wild_type/wildtype1_bacteriaGrow_cellsDie.tif");
		// Wild Type 2
		ImagePlus wt_2 = IJ.openImage(path + "wild_type/wildtype2_bacteriaGrow_cellsDie.tif");
		// Mutant
		ImagePlus mutant = IJ.openImage(path + "mutant/mutant_bacteriaGrowLess_cellsSurvive.tif");

		// 2. Process and segment the images
		// Wild Type 1
		ImagePlus[] segmentedImages_wt1 = processImage(wt_1);
		ImagePlus segmentedInf_wt1 = segmentedImages_wt1[0];
		segmentedInf_wt1.setTitle("Wild type 1 - Infection");
		ImagePlus segmentedMacro_wt1 = segmentedImages_wt1[1];
		segmentedMacro_wt1.setTitle("Wild type 1 - Macrophage");
		// Wild Type 2
		ImagePlus[] segmentedImages_wt2 = processImage(wt_2);
		ImagePlus segmentedInf_wt2 = segmentedImages_wt2[0];
		segmentedInf_wt2.setTitle("Wild type 2 - Infection");
		ImagePlus segmentedMacro_wt2 = segmentedImages_wt2[1];
		segmentedMacro_wt2.setTitle("Wild type 2 - Macrophage");
		// Mutant
		ImagePlus[] segmentedImages_mutant = processImage(mutant);
		ImagePlus segmentedInf_mutant = segmentedImages_mutant[0];
		segmentedInf_mutant.setTitle("Mutant - Infection");
		ImagePlus segmentedMacro_mutant = segmentedImages_mutant[1];
		segmentedMacro_mutant.setTitle("Mutant - Macrophage");

		// 3. Obtain values of interest
		// Wild Type 1
		double[][] inf_prog_wt1 = infection_progression(segmentedInf_wt1);
		double[][] inf_macro_wt1 = infection_by_macrophage(segmentedInf_wt1, segmentedMacro_wt1);
		double[] frac_wt1 = fraction_infected_cells(inf_macro_wt1[0], inf_macro_wt1[1]);
		double[] global_wt1 = global_proportion_cells_infection(segmentedInf_wt1, segmentedMacro_wt1);
		// Wild Type 2
		double[][] inf_prog_wt2 = infection_progression(segmentedInf_wt2);
		double[][] inf_macro_wt2 = infection_by_macrophage(segmentedInf_wt2, segmentedMacro_wt2);
		double[] frac_wt2 = fraction_infected_cells(inf_macro_wt2[0], inf_macro_wt2[1]);
		double[] global_wt2 = global_proportion_cells_infection(segmentedInf_wt2, segmentedMacro_wt2);
		// Mutant
		double[][] inf_prog_mutant = infection_progression(segmentedInf_mutant);
		double[][] inf_macro_mutant = infection_by_macrophage(segmentedInf_mutant, segmentedMacro_mutant);
		double[] frac_mutant = fraction_infected_cells(inf_macro_mutant[0], inf_macro_mutant[1]);
		double[] global_mutant = global_proportion_cells_infection(segmentedInf_mutant, segmentedMacro_mutant);


		// 4. Plot values
		plot_lines(frac_wt1, frac_wt2, frac_mutant, "Fraction of infected cells versus time", "Time", "Fraction of infected cells (0-1)");
		plot_lines(global_wt1, global_wt2, global_mutant, "Global pixels proportion of cells infection versus time", "Time", "Global pixels proportion of cells infection (0-1)");
		plot_dots(inf_prog_wt1[0], inf_prog_wt1[1], inf_prog_wt2[0], inf_prog_wt2[1], inf_prog_mutant[0], inf_prog_mutant[1], "Infections area versus time", "Time", "Infections area (pixels)");
		plot_dots(inf_macro_wt1[0], inf_macro_wt1[1], inf_macro_wt2[0], inf_macro_wt2[1], inf_macro_mutant[0], inf_macro_mutant[1], "Percentage area infected per cells versus time", "Time", "Percentage area infected per cells (%)");

		// Close ROI manager and result table
		RoiManager roiManager = RoiManager.getInstance();
		if (roiManager != null) {
			roiManager.close();
		}
		ResultsTable rt = ResultsTable.getResultsTable();
		if (rt != null) {
			IJ.selectWindow("Results");
			IJ.run("Close");
		}
	}


	/**
	 * Plots three lines on a graph representing wildtype 1, wild type 2, and mutant data.
	 *
	 * @param y_wt1    The y-values of the wildtype 1 data.
	 * @param y_wt2    The y-values of the wildtype 2 data.
	 * @param y_mutant The y-values of the mutant data.
	 * @param title    The title of the plot.
	 * @param xlabel   The label for the x-axis.
	 * @param ylabel   The label for the y-axis.
	 */
	private void plot_lines(double[] y_wt1, double[] y_wt2, double[] y_mutant, String title, String xlabel, String ylabel) {

		// Create time points values
		double[] timePoints = new double[y_wt1.length];

		for (int i = 0; i < timePoints.length; i++) {
			timePoints[i] = i + 1;
		}

		// Find the minimum and maximum y-values
		double minY = Math.min(Math.min(Arrays.stream(y_wt1).min().getAsDouble(), Arrays.stream(y_wt2).min().getAsDouble()), Arrays.stream(y_mutant).min().getAsDouble());
		double maxY = Math.max(Math.max(Arrays.stream(y_wt1).max().getAsDouble(), Arrays.stream(y_wt2).max().getAsDouble()), Arrays.stream(y_mutant).max().getAsDouble());

		// Create the plot
		Plot plot = new Plot(title, xlabel, ylabel);

		// Set the y-axis range to include all y-values in the window
		plot.setLimits(0, timePoints.length + 1, minY, maxY);

		// Add data points to the plot
		plot.setColor("red");
		plot.setLineWidth(2.0f);
		plot.addPoints(timePoints, y_wt1, Plot.LINE);
		plot.setColor("orange");
		plot.addPoints(timePoints, y_wt2, Plot.LINE);
		plot.setColor("blue");
		plot.addPoints(timePoints, y_mutant, Plot.LINE);

		// Create a legend
		plot.setColor(Color.RED);
		plot.addLabel(0.5, 0.25, "Wild Type 1");
		plot.setColor(Color.ORANGE);
		plot.addLabel(0.5, 0.2, "Wild Type 2");
		plot.setColor(Color.BLUE);
		plot.addLabel(0.5, 0.15, "Mutant");

		// Show the plot
		plot.show();
	}

	/**
	 * Plots three datasets from wildtype 1, wildtype 2, and mutant as dots on a graph.
	 * This method also add lines plots representing the mean of the data on each slice.
	 *
	 * @param x_wt1    The x-values of the wildtype 1 data.
	 * @param y_wt1    The y-values of the wildtype 1 data.
	 * @param x_wt2    The x-values of the wildtype 2 data.
	 * @param y_wt2    The y-values of the wildtype 2 data.
	 * @param x_mutant The x-values of the mutant data.
	 * @param y_mutant The y-values of the mutant data.
	 * @param title    The title of the plot.
	 * @param xlabel   The label for the x-axis.
	 * @param ylabel   The label for the y-axis.
	 */
	private void plot_dots(double[] x_wt1, double[] y_wt1, double[] x_wt2, double[] y_wt2, double[] x_mutant, double[] y_mutant, String title, String xlabel, String ylabel) {

		// Compute mean of x values
		double[] mean_wt1 = mean_values(x_wt1, y_wt1);
		double[] mean_wt2 = mean_values(x_wt2, y_wt2);
		double[] mean_mutant = mean_values(x_mutant, y_mutant);

		// Create time points values
		double[] timePoints = new double[mean_wt1.length];

		for (int i = 0; i < timePoints.length; i++) {
			timePoints[i] = i + 1;
		}


		// Create the plot
		Plot plot = new Plot(title, xlabel, ylabel);

		// Add data points to the plot
		plot.setColor("red");
		plot.addPoints(x_wt1, y_wt1, Plot.CIRCLE);
		plot.setColor("orange");
		plot.addPoints(x_wt2, y_wt2, Plot.CIRCLE);
		plot.setColor("blue");
		plot.addPoints(x_mutant, y_mutant, Plot.CIRCLE);

		// Add mean lines to the plot
		plot.setColor("red");
		plot.setLineWidth(2.0f);
		plot.addPoints(timePoints, mean_wt1, Plot.LINE);
		plot.setColor("orange");
		plot.addPoints(timePoints, mean_wt2, Plot.LINE);
		plot.setColor("blue");
		plot.addPoints(timePoints, mean_mutant, Plot.LINE);

		// Manually create a legend
		plot.setColor(Color.RED);
		plot.addLabel(0.5, 0.25, "Wild Type 1");
		plot.setColor(Color.ORANGE);
		plot.addLabel(0.5, 0.2, "Wild Type 2");
		plot.setColor(Color.BLUE);
		plot.addLabel(0.5, 0.15, "Mutant");

		// Show the plot
		plot.show();
	}

	/**
	 * The following method processes and segments the macrophage and infection channels.
	 *
	 * @param imp The ImagePlus image
	 * @return Two images corresponding to the infection and macrophage, processed and segmented.
	 */
	private ImagePlus[] processImage(ImagePlus imp) {

		// Slit the channels
		ImagePlus[] channels = ChannelSplitter.split(imp);
		ImagePlus infection = channels[0];
		ImagePlus macrophage = channels[1];

		// Set titles for the images
		infection.setTitle("Infection");
		macrophage.setTitle("Macrophage");
		infection.show();
		macrophage.show();


		// Step 1: Macrophage channel segmentation

		// Z projection on Max intensity
		macrophage = ZProjector.run(macrophage, "max all");
		// Subtract the background and filtering
		IJ.run(macrophage, "Subtract Background...", "rolling=50 sliding disable stack");
		IJ.run(macrophage, "Gaussian Blur...", "sigma=1 stack");
		IJ.run(macrophage, "Median...", "radius=5 stack");
		// Segmentation over variance
		IJ.run(macrophage, "Variance...", "radius=5 stack");
		// Gaussian filter and remove outliers
		IJ.run(macrophage, "Gaussian Blur...", "sigma=1 stack");
		IJ.run(macrophage, "Remove Outliers...", "radius=10 threshold=50 which=Bright stack");
		// Apply the segmentation
		IJ.run(macrophage, "Convert to Mask", "background=Dark calculate");
		// Fill holes and watershed
		IJ.run(macrophage, "Fill Holes", "stack");
		IJ.run(macrophage, "Watershed", "stack");

		// Show the segmented macrophage
		macrophage.setTitle("Segmented macrophage");
		macrophage.show();


		// Step 2: Infection channel segmentation

		// Z projection on Max intensity
		infection = ZProjector.run(infection, "max all");
		// Subtract background
		IJ.run(infection, "Subtract Background...", "rolling=50 stack");
		// Apply the segmentation
		IJ.run(infection, "Convert to Mask", "method=Otsu background=Dark calculate");

		// Show the segmented infection
		infection.setTitle("Segmented infection");
		infection.show();


		// Close previous images
		ImagePlus macro = WindowManager.getImage("Macrophage");
		if (macro != null) {
			macro.close();
		}
		ImagePlus inf = WindowManager.getImage("Infection");
		if (inf != null) {
			inf.close();
		}

		return new ImagePlus[]{infection, macrophage};
	}

	/**
	 * The following method calculates the global proportion of infected pixels over macrophages pixels for each time point.
	 * The algorithm uses the pixel's proportions of segmented infection and segmented cells obtained by calling
	 * the get_area method.
	 *
	 * @param segmented_inf   The ImagePlus object representing the segmented infection image.
	 * @param segmented_macro The ImagePlus object representing the segmented macrophages image.
	 * @return An array containing the global infection proportion for each time point.
	 */
	private double[] global_proportion_cells_infection(ImagePlus segmented_inf, ImagePlus segmented_macro) {

		// Compute the pixel's proportion of cells and infection at each time point
		List<Double> proportion_cell = get_proportion(segmented_macro);
		List<Double> proportion_inf = get_proportion(segmented_inf);

		// Calculate the pixels proportion for each time point
		List<Double> global_infection = new ArrayList<>();
		for (int i = 0; i < proportion_inf.size(); i++) {
			double infectionProportion = proportion_inf.get(i) / proportion_cell.get(i);
			global_infection.add(infectionProportion);
		}

		return global_infection.stream().mapToDouble(Double::doubleValue).toArray(); //Transform list into array
	}


	/**
	 * This method calculates the area of each segmented infection across all frames.
	 * The algorithm uses the segmented infection image to perform an analysis of particles and retrieve the ROIs area.
	 *
	 * @param segmented_inf The segmented image for measurements.
	 * @return A 2D double array containing slice numbers and infection area measurements.
	 */
	private double[][] infection_progression(ImagePlus segmented_inf) {

		// Apply the threshold and analyze particles
		IJ.run(segmented_inf, "Analyze Particles...", "exclude clear overlay add stack");

		// Ensure segmented_inf is the active image
		segmented_inf.show();

		// Initialize ROI Manager
		RoiManager rm = RoiManager.getRoiManager();
		Prefs.showAllSliceOnly = true;

		// Set Measurements for ROIs
		IJ.run("Set Measurements...", "area mean min area_fraction stack redirect=None decimal=3");

		// Perform measurements on the ROIs
		rm.runCommand(segmented_inf, "Measure");

		// Get the open ResultsTable
		ResultsTable rt = ResultsTable.getResultsTable();

		// Extract slice numbers and infection area measurements from the ResultsTable
		double[] slice = rt.getColumn("Slice");
		double[] area = rt.getColumn("Area");

		// Clear measurement results to avoid interference with other analyses
		IJ.run("Clear Results", "");

		// Deselect any selected ROIs in the ROI Manager
		rm.deselect();

		return new double[][]{slice, area};
	}

	/**
	 * The following method calculates the percentage of each macrophage area that is infected.
	 * The algorithm analyses the particles of the macrophages segmented image and retrieves the percentage of
	 * infected pixels inside the area of each ROI.
	 *
	 * @param segmented_inf   The segmented image of infections.
	 * @param segmented_macro The segmented image of macrophages.
	 * @return A 2D double array containing slice numbers and infected area percentages.
	 */
	private double[][] infection_by_macrophage(ImagePlus segmented_inf, ImagePlus segmented_macro) {

		// Analyze particles in the macrophage image
		IJ.run(segmented_macro, "Analyze Particles...", " exclude clear overlay add stack");

		// Initialize ROI Manager
		RoiManager rm = RoiManager.getRoiManager();
		Prefs.showAllSliceOnly = true;

		// Set Measurements for ROIs
		IJ.run("Set Measurements...", "area mean min area_fraction stack redirect=None decimal=3");

		// Run measurements on segmented_inf
		rm.runCommand(segmented_inf, "Measure");

		// Get the ResultsTable containing measurements
		ResultsTable rt = ResultsTable.getResultsTable();

		// Extract slice numbers and area percentages
		double[] slice = rt.getColumn("Slice");
		double[] perc_area = rt.getColumn("%Area");

		// Deselect ROIs
		rm.deselect();

		return new double[][]{slice, perc_area};

	}

	/**
	 * This method calculates the fraction of infected cells based on the given percentage of infected area for each cell
	 * and at each slice. A cell is considered infected if the infected area is greater than 0%
	 * Cells with infected area = 0% are not considered infected.
	 *
	 * @param slices The array containing slice numbers of each detected cell.
	 * @param ratio  The array containing %Area infection of each detected cell.
	 * @return An array containing the fraction of infected cells for each slice.
	 */
	private double[] fraction_infected_cells(double[] slices, double[] ratio) {

		// Map to store counts of total cells and infected cells > 1% area per slice
		Map<Integer, int[]> sliceCounts = new HashMap<>();

		for (int i = 0; i < slices.length; i++) {
			int slice = (int) slices[i];
			// Initialize counts for total cells and infected cells > 1% area
			sliceCounts.putIfAbsent(slice, new int[2]); // [total, %area>1]
			sliceCounts.get(slice)[0]++; // Increment total cell count
			if (ratio[i] > 0) sliceCounts.get(slice)[1]++; // Increment infected cell count (%area>1)
		}

		// Prepare the fraction of infected cells list
		List<Double> fractionList = new ArrayList<>();
		for (int i = 1; i <= sliceCounts.size(); i++) {
			int total = sliceCounts.getOrDefault(i, new int[]{0, 0})[0]; // Total cell count for the slice
			int infected = sliceCounts.getOrDefault(i, new int[]{0, 0})[1]; // Infected cell count for the slice
			fractionList.add(total > 0 ? (double) infected / total : 0); // Calculate fraction of infected cells
		}

		// Convert list to array for plotting
		double[] fraction = fractionList.stream().mapToDouble(Double::doubleValue).toArray();

		// Clear measurement results to avoid interference with other analyses
		IJ.run("Clear Results", "");

		return fraction;
	}

	/**
	 * This method calculates the proportion of black pixels in each slice of the given image.
	 * If applied on a segmented image, this method computes the pixel's proportion of segmented particles.
	 *
	 * @param imp The ImagePlus object representing the image.
	 * @return A list containing the proportion of non-black pixels in each slice of the image.
	 */
	private List<Double> get_proportion(ImagePlus imp) {

		// Get the number of slices in the hyperstack
		int numSlices = imp.getImageStackSize();

		// Initialize a list to store the proportions of black pixels for each slice
		List<Double> proportionsList = new ArrayList<>();

		// Iterate over each slice in the hyperstack
		for (int sliceIndex = 1; sliceIndex <= numSlices; sliceIndex++) {
			// Get the processor for the slice
			ImageProcessor sliceProcessor = imp.getImageStack().getProcessor(sliceIndex);

			// Get the statistics for the slice
			ImageStatistics stats = ImageStatistics.getStatistics(sliceProcessor, ImageStatistics.MIN_MAX, null);

			// Get the count of black pixels
			int blackPixelCount = stats.histogram[0];

			// Total number of pixels
			int totalPixels = sliceProcessor.getWidth() * sliceProcessor.getHeight();

			// Proportion of black pixels over white pixels
			double blackPixelProportion = (double) blackPixelCount / totalPixels;

			// Add the proportion to the list
			proportionsList.add(1 - blackPixelProportion);
		}

		return proportionsList;

	}

	/**
	 * This method calculates the mean measurement of each slice.
	 * As each slice don't necessary have the same number of measurements, this algorithm uses a list of slice
	 * in order to identify the slice associated to each measurement.
	 *
	 * @param slice The slice index of each measurement
	 * @param data  The array containing measurements
	 * @return A list containing the mean data of each slice.
	 */
	private double[] mean_values(double[] slice, double[] data) {

		// Get unique slice indices
		double[] unique_slices = java.util.Arrays.stream(slice).distinct().toArray();
		// Initialize array to store mean values
		double[] mean_by_slice = new double[unique_slices.length];
		for (int i = 0; i < unique_slices.length; i++) {
			// Current slice index
			double slice_nb = unique_slices[i];
			// Initialize sum and count of measurements for the current slice
			double sum = 0;
			int count = 0;
			for (int j = 0; j < slice.length; j++) {
				// Check if the measurement belongs to the current slice
				if (slice[j] == slice_nb) {
					// Add measurement to the sum
					sum += data[j];
					count++;
					// Increment count of measurements
				}
			}
			// Calculate the mean
			mean_by_slice[i] = sum / count;
		}
		return mean_by_slice;
	}


	/**
	 * This main function serves for development purposes.
	 * It allows you to run the plugin immediately out of
	 * your integrated development environment (IDE).
	 *
	 * @param args whatever, it's ignored
	 * @throws Exception
	 */
	public static void main(final String... args) throws Exception {
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();
	}
}

