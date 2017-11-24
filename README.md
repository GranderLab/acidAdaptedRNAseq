# acidAdaptedRNAseq

## Installation
1. Open R and install the R devtools package (if it is not already installed) with ```install.packages('devtools')```.
2. Download and install the acidAdaptedRNAseq package by entering ```install_github('GranderLab/acidAdaptedRNAseq')```.
3. Load the package with ```library(acidAdaptedRNAseq```.

## Usage
There is only one function included in the package; ```plotFigure``` which can be called with either the ```"go_analysis"```or ```"RNAseq"``` arguments. This function first replicates the operating system in which the analysis was performed (and thus can take some time to run) after which it renders the analysis file (.Rmd) and opens the resulting analysis in your web browser window. Here you can see the code and functions which have been used for each analysis and the resulting figures reproduced in real-time.

## Comments, Bugs, Suggestions
We are more than happy to help you with any questions you might have or hear your suggestions. Please submit these via the repositories issues section.

## Advanced usage
For the vast majority of users we imagine that the included ```plotFigure```function will be sufficient, although, all hand written functions used for analysis are included in the package. Any functions not included in external package dependencies are included in the R directory with the exception of some of the RNAseq analysis functions which are located in the inst/data-raw directory.