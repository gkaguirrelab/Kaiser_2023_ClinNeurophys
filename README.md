# Kaiser_2024_Psychophysiology
Data and analysis code pertaining to this publication:

A Kaiser, EB Haggerty, DP Garner, VY Bunya, GK Aguirre. (2024) A measure of the blink reflex to parametric variation of mechanical stimulation of the trigeminal nerve. Psychophysiology.

The lid position data for all subjects is present in the data directory. To reproduce the figures in the main body of the paper, run `primaryBlinkPCAnalysis.m` followed by `primaryBlinkPCAFigures.m`. The routine `tableBlinkRawVals.m` saves a table of measurement of aspects of the blink response that does not make use of the PCA analysis. Finally, the routine `preRegAnalysis.m` saves a table of results that reflect the original, pre-registered analysis strategy.

The paper features a supplementary materials section. The blink data presented in the supplement was collected and analyzed separately from the data in the main body of the paper. These data (and the analysis code) can be found in the repo [BLNK_2023_Expt](https://github.com/gkaguirrelab/BLNK_2023_Expt). There are measurements of the physical, acoustic properties of the stimulus. The measurements and plotting code are in the DropBox directory that GKA has for this paper, and thus not in this public repository.

