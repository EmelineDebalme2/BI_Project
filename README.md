<img src="icons/bii-banner.png" alt="image" width=auto height=40px>
<hr>

### **[Moodle BIO-410](https://moodle.epfl.ch/course/view.php?id=15721)**  <img src="icons/epfl-logo.png" width=auto height=40px alt="image" align="right">

><h3 style="font-weight:500; font-size:1.6em">Project INF : Infection Quantification Plugin</h3>
>
>


## Description
The goal of the project is to quantify the Mycobacterium tuberculosis infection in macrophages. Mycobacterium tuberculosis, i.e. Koch’s bacillus, responsible for tuberculosis, has been introduced in healthy macrophages and the cell-bacteria reaction has been imaged over one week long. Macrophages are immune cells which digest pathogens and foreign substances.

This plugin processes images to segment and analyze infection stages over time and outputs quantitative data regarding the infection progression and the fraction of infected cells.

## Installation

### Dependencies
Before using the Infection Quantification Plugin, ensure you have the following software installed:

- ImageJ
- Fiji, an ImageJ distribution
- Java Development Kit (JDK) 8 or higher

## Setup

### Clone the repository locally
1. Copy the repository URL from git
2. Open a git bash terminal on your computer, at the location where you would like to clone the repo.
3. Enter `` git clone `` and do a right click > paste to paste the repository URL.

<div align="center">
    <img src="icons/git-clone-command.png" width="80%">
</div>

### Open the project in intelliJ
1. Open your project
2. Under ``File -> Project structure``, select ``JDK 1.8`` and language level ``8``
3. Click on ``Apply`` and ``OK``


### Generate the jar file
Under the maven center, double-click on ``insall``. It generates a ``target`` repository inside which you'll find your ``.jar`` file. 
A ``Infection_quantification-1.0.0.jar`` file is already provided.

## Usage of the plugin

### Run the plugin
1. Click on ``Run ProjectCommand.java``
2. Go on the ImageJ window it opens
3. Under ``Plugins --> BII`` click on `Infection Quantification`
4. A dialog will appear prompting you to select the `data` path

or 
1. Open ImageJ
2. Open ``Infection_quantification-1.0.0.jar``
3. Restart ImageJ
3. Under ``Plugins --> BII`` click on `Infection Quantification`
4. A dialog will appear prompting you to select the `data` path


### Select `data` path 
The `data` folder should be composed of two sub-folders : 
- `mutant` with the mutant tif image named `mutant_bacteriaGrowLess_cellsSurvive.tif`
- `wild_type` with the two wild type tif images named `wildtype1_bacteriaGrow_cellsDie.tif`and `wildtype2_bacteriaGrow_cellsDie.tif`

### Image Processing
- The plugin will open, split the channels and display the images.
- It will then process these images to segment the infected cells and analyze the infection progression.

### Analyze Results

The plugin will generate several images and plots, including:
- Segmented macrophages and infections for each image (Wild type 1-2, Mutant)
- Plot of global proportion of cells infection versus time
- Plot of the percentage of area infected per cells versus time
- Plot of the infections area versus time
- Plot of the fraction of infected cells versus time 
