#power_webapp

##About

This will be a Shiny app that allows one to see at what fold change difference and base abundance they can detect differences (shown by a P-value heat map calculated by simulation), given a hypothetical number of samples and reads per sample.

This can be used for designing RNAseq experiments, metagenomic gene catalog experiments, and 16S rRNA seq experiments. For 16S rRNA seq experiments, there is the option to simulate data from a body site (based on data from the human microbiome project).

##Run instructions

To run the app, type the following command into R:
```
shiny::runApp()
```

##Dependencies
Required packages include shiny, ALDEx2, gplots, RColorBrewer


