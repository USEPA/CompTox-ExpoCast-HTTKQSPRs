# Supplemental Material

## Collaborative Evaluation of In Silico Predictions for High Throughput Toxicokinetics

### Abstract

High throughput toxicokinetic (HTTK) methods help to address chemical risk assessment data gaps but typically require chemical-specific in vitro measurements. As an alternative to measurement, in silico models for HTTK have been developed. We used in vivo data to evaluate the ]performance of a high throughput physiologically based toxicokinetic (HT-PBTK) model based on in silico parameters. There were three points of reference 1) in vitro measured TK data, 2) empirical TK parameters estimated from the in vivo data, and 3) randomized in vitro TK data (true values for the wrong chemical). Six quantitative structure property relationship (QSPR) models were used to estimate intrinsic hepatic clearance (Clint), fraction of chemical unbound in plasma (fup), and/or TK elimination half-life. A sensitivity analysis showed that HTTK-based predictions of early time points, including Cmax, are generally insensitive to Clint and fup parameters. Cmax could be predicted by in vitro data with RMSLE 0.84 and with QSPRs 0.83-1.1. Models had greater impact at later time points and the related AUC and steady-state concentration (Css). HT-PBTK based on in vitro measurements predicted AUC with RMSLE 1.3 and QSPRs ranged from 1.1-1.4. Based upon a separate analysis of 131 chemicals with in vitro measured fup and Clint for both human and rat, we estimate that interspecies differences inflate the RMSLE for Cmax by 0.43. A consensus QSPR generally outperformed other approaches. For novel compounds, a consensus QSPR approach may yield reasonable predictions of TK. 
### Data Analysis Vignettes

The following [R](https://cran.r-project.org/ "R") [vignettes](https://r-pkgs.org/vignettes.html "Vignettes") were written in [R Markdown](https://rmarkdown.rstudio.com/ "R Markdown") using [RStudio](https://posit.co/downloads/ "Download RStudio"):

1. [Loading Data Files](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-1-Data.Rmd "Loading Data Files").
2. [Level 1 Analysis](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-2-Level1Analysis.Rmd "Level 1 Analysis").
3. [Make CvT Predictions](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-3-MakeLevel2and3Preds.Rmd "Make CvT Predictions").
4. [Level 2 Analysis](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-4-Level2AnalyzePreds.Rmd "Level 2 Analysis").
5. [Make Level 2 Plots](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-5-Level2MakePlots.Rmd "Make Level 2 Plots").
6. [Level 3 Analysis](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-6-Level3Analysis.Rmd "Level 3 Analysis").
7. [Sensitivity Analysis](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-7-SensitivityAnalysis.Rmd "Sensitivity Analysis").
8. [Interspecies Analysis](https://github.com/USEPA/CompTox-ExpoCast-HTTKQSPRs/blob/main/TKQSPRs-8-InterspecesEvaluation.Rmd "Interspecies Analysis").

### Authors

#### Corresponding author 
John Wambaugh [wambaugh.john@epa.gov]

#### Full Author List

John Wambaugh<sup>1</sup> [Wambaugh.john@epa.gov]

Nisha Sipes<sup>1</sup> [Sipes.nisha@epa.gov]

Gilberto Padilla-Mercado<sup>1,2</sup>	[PadillaMercado.Gilberto@epa.gov]

Jon Arnot<sup>3</sup>	[jon@arnotresearch.com]

Linda Bertato<sup>4</sup> [l.bertato@uninsubria.it]

Trevor Brown<sup>3</sup>

Nicola Chirico<sup>4</sup> [nicola.chirico@uninsubria.it]

Christopher Cook<sup>1,2</sup> 

Daniel Dawson<sup>1,5</sup>	[ddawson@integral-corp.com]

Sarah Davidson-Fritz<sup>1</sup>	[DavidsonFritz.Sarah@epa.gov]

John DiBella<sup>6</sup>	[john.dibella@simulations-plus.com]

Stephen Ferguson<sup>7</sup>	[stephen.ferguson@nih.gov]

Rocky Goldsmith<sup>1,8</sup>	[rgoldsmith@congruencetx.com]

Chris Grulke<sup>1</sup>	

Richard Judson<sup>1</sup>	

Michael Lawless<sup>6</sup>	[mlawless@simulations-plus.com]

Kamel Mansouri<sup>9</sup>	[kamel.mansouri@nih.gov]

Grace Patlewicz<sup>1</sup>	[Patlewicz.Grace@epa.gov]

Ester Papa<sup>4</sup>	[ester.papa@uninsubria.it]

Prachi Pradeep<sup>1,10</sup>	

Alessandro Sangion<sup>3</sup>	[alessandro.sangion@mail.utoronto.ca]

Risa Sayre<sup>1</sup>	[sayre.risa@epa.gov]

Russell Thomas<sup>1</sup>	[Thomas.Russell@epa.gov]

Rogelio Tornero-Velez<sup>1</sup>	[tornero-velez.rogelio@epa.gov]

Barbara Wetmore<sup>1</sup>	[wetmore.barbara@epa.gov]

Michael Devito<sup>1</sup>	[Devito.Michael@epa.gov]

1.	Center for Computational Toxicology and Exposure, Office of Research and Development, United States Environmental Protection Agency, Research Triangle Park, North Carolina, USA 28311
2.	Oak Ridge Institute for Science and Education, Oak Ridge, Tennessee, USA 38331
3.	ARC Arnot Research and Consulting Inc., Toronto, Ottawa, Canada
4.	Department of Theoretical and Applied Sciences, University of Insubria, Varese, Italy
5.	Integral Consulting, Seattle, Washington, USA 98104
6.	Simulations Plus, Inc., Lancaster, California, USA 93534
7.	Division of Translational Toxicology (DTT), National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina, USA 27709
8.	Congruence Therapeutics, Chapel Hill, North Carolina, USA 27517
9.	NTP Interagency Center for the Evaluation of Alternative Toxicological Methods, National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina, USA 28317
10.	German Federal Institute for Risk Assessment (BfR), Berlin-Jungfernheide, Germany

### License

License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>