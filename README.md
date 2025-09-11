Hi there, 

This is a pipeline to identify cell populations based on your cellomics data.

 1) Image your well plates as usual using the CX7.
 2) Upload the images to Columbus and use the attached analysis example to export the single cell data.
    Extract any information that could help in identifying neurons (e.g. Hoechst intensity, nucleus size, nucleus roundness, MAP2 intensity), 
    the transfection status (e.g. GFP intensity) and the intensity of your protein of interest (e.g. Drebrin).
 3) Put all the single cell data files into a folder.
 4) Open the R script and adapt it to your variables and number of clusters you can identify.
 5) Perform your statistical test of choice - perferably a multilevel analysis as the observations are not independent.

Have fun and reach out if you need help, 
Panthea

