Natalia Bailey, George Kapetanios, and M. Hashem Pesaran (2021),
"Measurement of Factor Strength: Theory and Practice",
Journal of Applied Econometrics, forthcoming.

Monte Carlo simulations (procedures and results)
(Section 5 of BKP 2021 paper and Supplementary Appendix E)

All Gauss procedures corresponding to the experiments of Section 5 of BKP 2021 main paper
can be found in Procedures\Simulations. These are:

File: "mc_FS1A.prg"
File: "mc_FS1A_Power.prg"
File: "mc_FS1B.prg"
File: "mc_FS1B_Power.prg"
File: "mc_FS2A.prg"
File: "mc_FS2A_Power.prg"
File: "mc_FS2B.prg"
File: "mc_FS2B_Power.prg"
File: "mc_FS3Aa.prg"
File: "mc_FS3Aa_Power.prg"
File: "mc_FS3Ab.prg"
File: "mc_FS3Ab_Power.prg"
File: "mc_FS3B_PC.prg"
File: "mc_FS3Ba_CSA.prg"
File: "mc_FS3Bb_CSA.prg"
File: "mc_FS3Bc_SeqCSA.prg"
File: "mc_FS4a.prg"
File: "mc_FS4b.prg"

The names of these files correspond to the names of the Experiments in Section 5 of BKP 2021 main paper.
The description of each experiment and the corresponding Tables/Figures in the paper are provided at 
the top of each .prg file.

The content of all results tables and figures to the experiments of Section 5 of BKP 2021
shown in the main paper and Supplementary Appendix E can be found in Results\Simulations. 
These are:

File: FS_MC_tables.xlsx
File: power_matlab.xlsx


Empirical Applications (procedures and results)
(Section 6 of BKP 2021 paper and Supplementary Appendix D)

The Gauss procedures corresponding to the financial and macroeconomic applications
of Section 6 of BKP 2021 main paper can be found in Procedures\Empirics. These are:

File: "aFS_finance.prg"
File: "aFS_macro.prg"

Data Description

Subsection 6.1: Identifying risk factors in asset pricing models

The S&P500 data (345 xls files for returns plus 3 files that include
market, risk factor and interest rate data) can be found in Procedures\Empirics.
These data read into the Gauss file: "aFS_finance.prg"

Subsection 6.2: Stength of common macroeconomic shocks

The macroeconomic data (FRED_QSWf.xls) can be found in Procedures\Empirics.
The raw data, which include both high-level economic and financial aggregates 
as well as disaggregated components, can be found on the Federal Reserve Bank
of St Louis website at: https://research.stlouisfed.org/econ/mccracken/static.html.
The data in file: "FRED_QSWf.xls" read into the Gauss file: "aFS_macro.prg"

Output

Output from "aFS_finance.prg" for the financial empirical application of Sub-section 6.1
of BKP 2021 paper can be found in Results\Empirics. All results and further computations
can be found in file: "FS_empirics_tables_figures.xlsx". Each sheet is named after the 
table or figure in the main BKP 2021 paper or Supplementary Appendix D that the results
correspond to.

Output from "aFS_macro.prg" corresponds directly to the content of Table 5 of BKP2021
main paper.


Please address any questions to:

Natalia Bailey
Department of Econometrics and Business Statistics,
Monash University, Melbourne, Australia
Email: natalia.bailey@monash.edu