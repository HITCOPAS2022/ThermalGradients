# ThermalGradient

This workflow starts by extracting daily sea surface temperature (SST) rasters from the ERDDAP server through the R package rerddapXtracto (Mendelssohn, 2020). Based on these data we use the Contextual Median Filter algorithm for calculating thermal gradients (TG) through the R package grec (Lau-Medrano, 2020). Original TG maps are rescaled to 10 × 10 km grid size raster, and values averaged over these new grid-cells. As an example, we conduct these analyses for the Northern Chilean Patagonia using SST data from January 2022. 

This is part of the research project "Coupling oceanographic and habitat models to assess abundance, distribution, and risk for baleen whales in Chile: Developing new tools for management” undertaken by researchers from Centro de Investigación Oceanográfica en el Pacífico Sur Oriental (COPAS Coastal), Universidad Austral de Chile, Universidad de Valparaíso, Universidad de Concepción, Instituto de Fomento Pesquero (IFOP), Centro de Estudios Avanzados en Zonas Áridas (CEAZA), Oregon State University, and Instituto Aqualie. Funded by COPAS Coastal HIT projects 2022.
 
![map](https://user-images.githubusercontent.com/14979334/228365844-82e7af30-9ced-4e83-87c7-d296df407d19.gif)
