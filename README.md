# Supplemental Material

## Collaborative Evaluation of In silico Predictions for High Throughput Toxicokinetics

### Abstract

To assess public health risks posed by chemicals we need chemical-specific toxicokinetic (TK) data to understand chemical absorption, distribution, metabolism, and elimination by the body. High throughput TK (HTTK) methods have the potential to address data gaps. This collaborative trial uses a database of in vivo measured toxicokinetic data to evaluate in silico approaches for HTTK. Six different sets of quantitative structure-property-relationship (QSPR) tools were evaluated. The predicted parameter values were used within a high throughput physiologically based TK (PBTK) model to predict in vivo measured plasma concentrations for 92 chemicals, mostly in rats. Root mean squared log10 error (RMSLE) was calculated for toxicologically-relevant dosimetry statistics AUC (time-integrated area under the curve) and Cmax (peak concentration). Early time points, including Cmax, are driven by physicochemical properties and are insensitive to in vitro measured HTTK data. Cmax could be predicted with RMSLE 0.8-0.9 in contrast to 0.7 for in vivo fits. For AUC in vivo parameters gave RMSLE 0.5. Chemical-specific in vitro data had RMSLE 1.3 and QSPRs ranged from 1.3-1.4. In vivo parameters for the full concentration time-course data had RMSLE of 0.62. In vitro HTTK predicted gave RMSLE of 1.18. Greater discrimination between QSPRs was observed at later time points, which impacts AUC and is driven by estimated metabolism and elimination. QSPRs predict plasma binding well (RMSLE 0.03 – 0.07) but have difficulty predicting metabolic clearance (RMSLE 0.37 – 1.28). A consensus prediction using the maximum clearance predicted across all QSPRs predicted AUC with RMSLE 1.1 and the full time-course with RMSLE 1.1. The consensus predictions outperformed the in vitro measured data for the evaluation chemicals. For novel compounds a consensus QSPR approach may yield reasonable predictions of TK. This evaluation characterizes the accuracy of HTTK approaches for new chemicals based on both in vitro measurement and structure-based in silico predictions.

### Data Analysis Vignettes

The following [R](https://cran.r-project.org/ "R") [vignettes](https://r-pkgs.org/vignettes.html "Vignettes") were written in [R Markdown](https://rmarkdown.rstudio.com/ "R Markdown") using [RStudio](https://posit.co/downloads/ "Download RStudio"):

1. [Loading Data Files](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-1-Data.Rmd "Loading Data Files").
2. [Level 1 Analysis](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-2-Level1Analysis.Rmd "Level 1 Analysis").
3. [Make CvT Predictions](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-3-Level2MakePreds.Rmd "Make CvT Predictions").
4. [Level 2 Analysis](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-4-Level2AnalyzePreds.Rmd "Level 2 Analysis").
5. [Make Level 2 Plots](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-5-Level2MakePlots.Rmd "Make Level 2 Plots").
6. [Level 3 Analysis](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-6-Level3Analysis.Rmd "Level 3 Analysis").

### Authors

#### Principal Investigator 
John Wambaugh [wambaugh.john@epa.gov]

#### Full Author List

John Wambaugh<sup>1</sup> [Wambaugh.john@epa.gov]

Nisha Sipes<sup>1</sup> [Sipes.nisha@epa.gov]

Kamel Mansouri<sup>2,3</sup>	[kamel.mansouri@nih.gov]

Jon Arnot<sup>4</sup>	[jon@arnotresearch.com]

Trevor Brown<sup>4</sup>	

Christopher Cook<sup>1</sup>	

Daniel Dawson<sup>1</sup>	[ddawson@integral-corp.com]

Sarah Davidson-Fritz<sup>1</sup>	[DavidsonFritz.Sarah@epa.gov]

John DiBella<sup>5</sup>	[john.dibella@simulations-plus.com]

Stephen Ferguson<sup>6</sup>	[stephen.ferguson@nih.gov]

Rocky Goldsmith<sup>1</sup>	[rgoldsmith@congruencetx.com]

Chris Grulke<sup>1</sup>	

Richard Judson<sup>1</sup>	[Judson.Richard@epa.gov]

Michael Lawless<sup>5</sup>	[mlawless@simulations-plus.com]

Gilberto Padilla-Mercado<sup>1</sup>	[PadillaMercado.Gilberto@epa.gov]

Grace Patlewicz<sup>1</sup>	[Patlewicz.Grace@epa.gov]

Ester Papa<sup>7</sup>	[ester.papa@uninsubria.it]

Prachi Pradeep<sup>8</sup>	

Alessandro Sangion<sup>2</sup>	[alessandro.sangion@mail.utoronto.ca]

Risa Sayre<sup>1</sup>	[sayre.risa@epa.gov]

Russell Thomas<sup>1</sup>	[Thomas.Russell@epa.gov]

Rogelio Tornero-Velez<sup>1</sup>	[tornero-velez.rogelio@epa.gov]

Barbara Wetmore<sup>1</sup>	[wetmore.barbara@epa.gov]

Michael Devito<sup>1</sup>	[Devito.Michael@epa.gov]

1.	U.S. Environmental Protection Agency, Office of Research and Development, Center for Computational Toxicology and Exposure, Research Triangle Park, NC 27711, USA
2. Innotiv
3. NICEATM
4. Arnot Research & Consulting
5. Simulations Plus
6. NIEHS
7. University of Insubria, Varese
8. German Federal Institute for Risk Assessment (BfR)

### License

License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>