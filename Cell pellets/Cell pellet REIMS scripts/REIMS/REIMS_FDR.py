'''FDR for the REIMS data. This code applies FDR correction using the Benjamini-Hochberg Procedure.'''

'''___Authorship___'''
__author__ = "Ariadna Gonzalez"
__created__ = "June 2023"
__version__ = "v1"
__maintainer__ = "Ariadna Gonzalez"
__email__ = "ariadna.gonzalez@npl.co.uk"

'''___Third-Party Modules___'''
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as sm

# Define inputs and outputs path
path_data = '...'
path_outputs = '...'


# Read excel file as a pandas DataFrame
df_volcano = pd.read_excel(path_data + "PariwiseTTestMaleVsFemaleWT.xlsx")

#Define column of interest and apply FDR
fdr_p = df_volcano["p-value Male to Female"]
fdr_p = np.asarray(fdr_p)
fdr, p_corrected,alpha1,alpha2 = sm.multipletests(fdr_p, 0.05, method='fdr_bh', is_sorted=False)

#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr = df_volcano.copy()
df_volcano_fdr.insert(2, 'fdr_bh Male to Female', fdr)
df_volcano_fdr.insert(3, "adjusted p-value Male to Female", p_corrected)

#Save new excel with adjusted p-values in outputs path
df_volcano_fdr.to_excel(path_outputs + "PariwiseTTestMaleVsFemaleWT_adjusted.xlsx", index=False)




# Read the excel file as a pandas DataFrame
df_volcano = pd.read_excel(path_data + "epithelialREIMScellTtestLog2geneticVariants.xls")

#Define column of interest and apply FDR
fdr_p = df_volcano["p-value APC KRAS PTEN to APC"]
fdr_p = np.asarray(fdr_p)
fdr, p_corrected,alpha1,alpha2 = sm.multipletests(fdr_p, 0.05, method='fdr_bh', is_sorted=False)

#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr = df_volcano.copy()
df_volcano_fdr.insert(2, 'fdr_bh APC KRAS PTEN to APC', fdr)
df_volcano_fdr.insert(3, "adjusted p-value APC KRAS PTEN to APC", p_corrected)

#Define column of interest and apply FDR
fdr_p2 = df_volcano["p-value APC KRAS PTEN to APC KRAS"]
fdr_p2 = np.asarray(fdr_p2)
fdr2, p_corrected2,alpha1,alpha2 = sm.multipletests(fdr_p2,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(7, 'fdr_bh APC KRAS PTEN to APC KRAS', fdr2)
df_volcano_fdr.insert(8, "adjusted p-value APC KRAS PTEN to APC KRAS", p_corrected2)

#Define column of interest and apply FDR
fdr_p3 = df_volcano["p-value APC KRAS PTEN to WT"]
fdr_p3 = np.asarray(fdr_p3)
fdr3, p_corrected3,alpha1,alpha2 = sm.multipletests(fdr_p3,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(12, 'fdr_bh APC KRAS PTEN to WT', fdr3)
df_volcano_fdr.insert(13, "adjusted p-value APC KRAS PTEN to WT", p_corrected3)

#Define column of interest and apply FDR
fdr_p4 = df_volcano["p-value APC KRAS PTEN to KRAS"]
fdr_p4 = np.asarray(fdr_p4)
fdr4, p_corrected4,alpha1,alpha2 = sm.multipletests(fdr_p4,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(17, 'fdr_bh APC KRAS PTEN to KRAS', fdr4)
df_volcano_fdr.insert(18, "adjusted p-value APC KRAS PTEN to KRAS", p_corrected4)

#Define column of interest and apply FDR
fdr_p5 = df_volcano["p-value APC to APC KRAS"]
fdr_p5 = np.asarray(fdr_p5)
fdr5, p_corrected5,alpha1,alpha2 = sm.multipletests(fdr_p5,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(22, 'fdr_bh APC to APC KRAS', fdr5)
df_volcano_fdr.insert(23, "adjusted p-value APC to APC KRAS", p_corrected5)

#Define column of interest and apply FDR
fdr_p6 = df_volcano["p-value APC to WT"]
fdr_p6 = np.asarray(fdr_p6)
fdr6, p_corrected6,alpha1,alpha2 = sm.multipletests(fdr_p6,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(27, 'fdr_bh APC to WT', fdr6)
df_volcano_fdr.insert(28, " adjusted p-value APC to WT", p_corrected6)

#Define column of interest and apply FDR
fdr_p7 = df_volcano["p-value APC to KRAS"]
fdr_p7 = np.asarray(fdr_p7)
fdr7, p_corrected7,alpha1,alpha2 = sm.multipletests(fdr_p7,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(32, 'fdr_bh APC to KRAS', fdr7)
df_volcano_fdr.insert(33, "adjusted p-value APC to KRAS", p_corrected7)

#Define column of interest and apply FDR
fdr_p8 = df_volcano["p-value APC KRAS to WT"]
fdr_p8 = np.asarray(fdr_p8)
fdr8, p_corrected8,alpha1,alpha2 = sm.multipletests(fdr_p8,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(37, 'fdr_bh APC KRAS to WT', fdr8)
df_volcano_fdr.insert(38, "adjusted p-value APC KRAS to WT", p_corrected8)

#Define column of interest and apply FDR
fdr_p9 = df_volcano["p-value APC KRAS to KRAS"]
fdr_p9 = np.asarray(fdr_p9)
fdr9, p_corrected9,alpha1,alpha2 = sm.multipletests(fdr_p9,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(42, 'fdr_bh APC KRAS to KRAS', fdr9)
df_volcano_fdr.insert(43, "adjusted p-value APC KRAS to KRAS", p_corrected9)

#Define column of interest and apply FDR
fdr_p10 = df_volcano["p-value WT to KRAS"]
fdr_p10 = np.asarray(fdr_p10)
#Insert columns with fdr and adjusted p-values in DataFrame
fdr10, p_corrected10,alpha1,alpha2 = sm.multipletests(fdr_p10,0.05, method='fdr_bh', is_sorted=True)
df_volcano_fdr.insert(47, 'fdr_bh WT to KRAS', fdr10)
df_volcano_fdr.insert(48, "adjusted p-value WT to KRAS", p_corrected10)

#Save new excel with adjusted p-values in outputs path
df_volcano_fdr.to_excel(path_outputs + "epithelialREIMScellTtestLog2geneticVariants_adjusted.xlsx", index=False)




# Read the excel file as a pandas DataFrame
df_volcano = pd.read_excel(path_data + "maleVsFemaleTTestLog2Genotypes.xlsx")

#Define column of interest and apply FDR
fdr_p = df_volcano["p-value APC, Kras, Pten male to APC, Kras, Pten female"]
fdr_p = np.asarray(fdr_p)
fdr, p_corrected,alpha1,alpha2 = sm.multipletests(fdr_p, 0.05, method='fdr_bh', is_sorted=False)

#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr = df_volcano.copy()
df_volcano_fdr.insert(2, 'fdr_bh APC, Kras, Pten male to APC, Kras, Pten female', fdr)
df_volcano_fdr.insert(3, "adjusted p-value APC, Kras, Pten male to APC, Kras, Pten female", p_corrected)

#Define column of interest and apply FDR
fdr_p2 = df_volcano["p-value APC, Kras, male to APC, Kras female"]
fdr_p2 = np.asarray(fdr_p2)
fdr2, p_corrected2,alpha1,alpha2 = sm.multipletests(fdr_p2,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(7, 'fdr_bh APC, Kras, male to APC, Kras female', fdr2)
df_volcano_fdr.insert(8, "adjusted p-value APC, Kras, male to APC, Kras female", p_corrected2)

#Define column of interest and apply FDR
fdr_p3 = df_volcano["p-value APC male to APC female"]
fdr_p3 = np.asarray(fdr_p3)
fdr3, p_corrected3,alpha1,alpha2 = sm.multipletests(fdr_p3,0.05, method='fdr_bh', is_sorted=True)
#Insert columns with fdr and adjusted p-values in DataFrame
df_volcano_fdr.insert(12, 'fdr_bh APC male to APC female', fdr3)
df_volcano_fdr.insert(13, "adjusted p-value APC male to APC female", p_corrected3)

#Save new excel with adjusted p-values in outputs path
df_volcano_fdr.to_excel(path_outputs + "maleVsFemaleTTestLog2Genotypes_adjusted.xlsx", index=False)
