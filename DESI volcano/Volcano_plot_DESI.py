'''Volvano plot DESI analysis.
This code performs the t-test on the APC vs APC/KRAS DESI samples. It applies FDR correction with the
Benjamini-Hochberg Procedure.'''

'''___Authorship___'''
__author__ = "Ariadna Gonzalez"
__created__ = "June 2023"
__version__ = "v2"
__maintainer__ = "Ariadna Gonzalez"
__email__ = "ariadna.gonzalez@npl.co.uk"

'''___Third-Party Modules___'''
import numpy as np
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, normaltest, norm
import statsmodels.stats.multitest as sm

path_data = "..."
path_outputs= "..."

# Data for each of tissue-masks is obtained from the datacube () and the masks () provided
mask_A_apckras = np.genfromtxt(path_data + "mask_A_apckras.csv",delimiter=",")
mask_D_apckras = np.genfromtxt(path_data+"mask_D_apckras.csv",delimiter=",")
mask_C_apc = np.genfromtxt(path_data+"mask_C_apc.csv",delimiter=",")
mask_G_apc = np.genfromtxt(path_data+"mask_G_apc.csv",delimiter=",")

#Concatenation of similar tumour types
mask_D_apckras = mask_D_apckras[:,1:]
mask_G_apc = mask_G_apc[:,1:]
apckras = np.concatenate((mask_A_apckras,mask_D_apckras), axis = 1)
apc = np.concatenate((mask_C_apc,mask_G_apc), axis = 1)

#Two tailed t-test for apc and apckras samples.
print("T-test running...")
ttest_table = np.zeros((len(apckras[:,0]),9))
for i in range(len(apckras[:,0])):
    apckras_intensities = apckras[i,1::]
    apckras_intensities = apckras_intensities[apckras_intensities != 0]
    mean_apckras_intensities = np.mean(apckras_intensities)

    apc_intensities = apc[i,1::]
    apc_intensities = apc_intensities[apc_intensities != 0]
    mean_apc_intensities = np.mean(apc_intensities)

    t=np.zeros(150)
    p=np.zeros(150)
    for jj in range(0, 150, 1):
        pixels_apckras_intensities = np.random.choice(apckras_intensities, 200, replace=True)
        pixels_apc_intensities = np.random.choice(apc_intensities, 200, replace=True)
        t,p = ttest_ind(pixels_apc_intensities,pixels_apckras_intensities,equal_var=False)
    t_mean = np.mean(t)
    p_mean = np.mean(p)
    t_std = np.std(t)
    p_std = np.std(p)

    ttest_table[i] = ((apckras[i,0],
                      mean_apc_intensities,
                      mean_apckras_intensities,
                      t_mean, p_mean, t_std,p_std,
                      -np.log10(p_mean),
                      np.log2(np.mean(apckras_intensities)/np.mean(apc_intensities))))
    
    #Save results in DataFrame
    df_volcano = pd.DataFrame(ttest_table,columns = ["m/z","intensity apc","intensity apckras",
                                                     "t_mean","p_mean","t_std",
                                                     "p_std","-log10(p)","log2(FC)"])

df_volcano['p_values'] = 10 ** -(df_volcano["-log10(p)"])
df_volcano_fdr = df_volcano.copy()

#Correcting FDR with Benjamini-Hochberg procedure
fdr_p = df_volcano["p_values"]
fdr_p = np.asarray(fdr_p)
fdr, p_corrected, alpha1, alpha2 = sm.multipletests(fdr_p, 0.05, method='fdr_bh', is_sorted=False)

#Insert outputs of FDR
df_volcano_fdr.insert(5, 'fdr_bh APC to APC/KRAS', fdr)
df_volcano_fdr.insert(6, "adjusted p-value APC to APC/KRAS", p_corrected)
df_volcano_fdr.insert(7, "-log10(adjusted p-value)", -np.log10(p_corrected))
df_volcano_fdr['sign adjusted'] = 'not significant'

#Significance criteria
df_volcano_fdr.loc[(df_volcano_fdr['-log10(adjusted p-value)'] > -np.log10(0.05)) & (
        df_volcano_fdr['log2(FC)'] > np.log2(1.5)), 'sign adjusted'] = 'up regulated in apc-kras'
df_volcano_fdr.loc[(df_volcano_fdr['-log10(adjusted p-value)'] > -np.log10(0.05)) & (
        df_volcano_fdr['log2(FC)'] < -np.log2(1.5)), 'sign adjusted'] = 'down regulated in apc-kras'

df_volcano_fdr.to_excel(path_outputs + "APC_APCKRAS_ttest_fdr_q0.05.xlsx", index=False)
print("!Results saved in: " + path_outputs)

#Creating volcano plot
print('Creating volcano plot...')
xaxis_lim = np.absolute(df_volcano_fdr['log2(FC)']).max(axis=0)*1.1
fig = px.scatter(df_volcano_fdr,x='log2(FC)',y='-log10(adjusted p-value)',color='sign adjusted',
                 hover_data=['m/z'],color_discrete_map={'not significant':'grey',
        'up regulated in apc-kras':'blue','down regulated in apc-kras':'red'},
                 labels={"sign adjusted p-value":"",
                         "not significant":"Not Significant",
                         "up regulated in apc-kras":"Up regulated in APC/KRAS",
                        'down regulated in apc-kras':'Down regulated in APC/KRAS'})

fig.update_xaxes(range=[-xaxis_lim,xaxis_lim])
fig.update_layout(font=dict(size = 18))
fig.add_hline(y=-np.log10(0.05), line_width=1, line_dash="dash", line_color="grey")
fig.add_vline(x=-np.log2(1.5), line_width=1, line_dash="dash", line_color="grey")
fig.add_vline(x=np.log2(1.5), line_width=1, line_dash="dash", line_color="grey")
fig.show()
fig.write_html(path_outputs + "volcano_APC_APCKRAS_ttest_fdr.html")
fig.write_image(path_outputs + "volcano_APC_APCKRAS_ttest_fdr.png")
print("!Figure saved in: " + path_outputs)

print("!Finished")