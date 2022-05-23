import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics
import warnings
warnings.simplefilter('ignore')

PATH_TO_XCR_TABLE = '/home/akosenkov/key/Bases/XCR_exist_chain.txt'

def calc_fraction_groups(samples, chain):
    """
    
    Calc Extra_Hyperexpanded, Hyperexpanded, Extra_Large, Large, Medium, Small
    Extra_Small, Rare, Extra_Rare fraction groups for required samples.
    :param samples: list, file names
    :param chain: str, required chain, can contain: alpha, beta, gamma, delta,
                  heavy, kappa, lambda, light, IgM, IgG
    :return : f_extra_rare, f_rare, f_extra_small, f_small, f_med, f_large,
              f_extra_large, f_hyp, f_extra_hyp
    
    """

    fractions_list = []
    # data collection for each sample
    for samp in samples:
        freqs = {
            'Extra_Rare': 0,
            'Rare': 0,
            'Extra_Small': 0,
            'Small': 0,
            'Medium': 0,
            'Large': 0,
            'Extra_Large': 0,
            'Hyperexpanded': 0,
            'Extra_Hyperexpanded': 0
        }

        # loading table and chain restriction
        if type(samp) == str:
            data = pd.read_csv(samp, sep='\t')
        else:
            data = samp

        if chain == 'light':
            data = data[list(
                np.array(data.chain == 'lambda')
                | np.array(data.chain == 'kappa'))].reset_index(drop=True)
        elif 'Ig' in chain:
            data = data[np.array(data.isotype == chain)].reset_index(drop=True)
        else:
            data = data[data.chain == chain].reset_index(drop=True)

        # fraction calculation if chain is 'light', 'IgG', 'IgM'
        if chain == 'light' or 'Ig' in chain:
            s = data.clonotype_count.sum()
            data.fraction = data.clonotype_count.apply(lambda x: x / s)

        # data separating by fraction
        freqs['Extra_Rare'] += data[(data.fraction <= 0.00000001)].fraction.sum()
        freqs['Rare'] += data[(data.fraction <= 0.0000001)
                           & (data.fraction > 0.00000001)].fraction.sum()
        freqs['Extra_Small'] += data[(data.fraction <= 0.000001)
                                  &
                                  (data.fraction > 0.0000001)].fraction.sum()
        freqs['Small'] += data[(data.fraction <= 0.00001)
                            & (data.fraction > 0.000001)].fraction.sum()
        freqs['Medium'] += data[(data.fraction <= 0.0001)
                             & (data.fraction > 0.00001)].fraction.sum()
        freqs['Large'] += data[(data.fraction <= 0.001)
                            & (data.fraction > 0.0001)].fraction.sum()
        freqs['Extra_Large'] += data[(data.fraction <= 0.01)
                                  & (data.fraction > 0.001)].fraction.sum()
        freqs['Hyperexpanded'] += data[(data.fraction <= 0.1)
                                    & (data.fraction > 0.01)].fraction.sum()
        freqs['Extra_Hyperexpanded'] += data[(data.fraction > 0.1)].fraction.sum()

        fractions_list.append(freqs)

    f_extra_rare = [samp['Extra_Rare'] for samp in fractions_list]
    f_rare = [samp['Rare'] for samp in fractions_list]
    f_extra_small = [samp['Extra_Small'] for samp in fractions_list]
    f_small = [samp['Small'] for samp in fractions_list]
    f_med = [samp['Medium'] for samp in fractions_list]
    f_large = [samp['Large'] for samp in fractions_list]
    f_extra_large = [samp['Extra_Large'] for samp in fractions_list]
    f_hyp = [samp['Hyperexpanded'] for samp in fractions_list]
    f_extra_hyp = [samp['Extra_Hyperexpanded'] for samp in fractions_list]
    return f_extra_rare, f_rare, f_extra_small, f_small, f_med, f_large, f_extra_large, f_hyp, f_extra_hyp

def fraction_healthy_cohort(chain, material_type, healthy_base, seq_type):
    '''
    
    Calc Hyperexpanded, Large, Medium, Small, Rare groups for healthy cohort (Target).
    :param chain: str, required chain, can contain: alpha, beta, gamma, delta, heavy,
                  kappa, lambda, light, IgM, IgG
    :param material_type: str, tissue type: PBMC or whatever you want (for 'Target_Lab'
                          if not PBMC, then T_cells or B_cells samples will be used).
    :param healthy_base: str, base for healthy cohort.
    :return : h_f_extra_rare, h_f_rare, h_f_extra_small, h_f_small, h_f_med, h_f_large,
              h_f_extra_large, h_f_hyp, h_f_extra_hyp
    
    '''
    # loading path table
    xcr_table = pd.read_csv(PATH_TO_XCR_TABLE, sep='\t')
    if chain in ['alpha', 'beta', 'gamma', 'delta']:
        file = 'TCR_qc_pass.txt'
        ch = chain
        if healthy_base == 'Target_Lab':
            if material_type == 'PBMC':
                t_m = 'PBMC'
            else:
                t_m = 'T_cells'
        else:
            t_m = material_type
    elif chain in ['heavy', 'kappa', 'lambda', 'light'] or 'Ig' in chain:
        file = 'BCR_qc_pass.txt'
        if 'Ig' in chain:
            ch = 'heavy'
        else:
            ch = chain
        if healthy_base == 'Target_Lab':
            if material_type == 'PBMC':
                t_m = 'PBMC'
            else:
                t_m = 'B_cells'
        else:
            t_m = material_type
    else:
        print('ERROR: WRONG CHAIN')

    # path table restriction by 'healthy_base', 'chain', 'material_type'
    if chain=='light':
        xcr_table = xcr_table[(xcr_table.Base == healthy_base)
                              & (xcr_table.Cell_type == t_m.replace(' ', '_'))
                              & (xcr_table.Type == seq_type)
                              & (xcr_table.Diagnosis == 'Healthy')
                              &
                              (xcr_table.Stimulation == 'Without_stimulation')]
        xcr_table = xcr_table[xcr_table['lambda'] | xcr_table.kappa]
    else:
        xcr_table = xcr_table[(xcr_table.Base == healthy_base)
                              & (xcr_table.Cell_type == t_m.replace(' ', '_'))
                              & (xcr_table.Type == seq_type)
                              & (xcr_table.Diagnosis == 'Healthy')
                              &
                              (xcr_table.Stimulation == 'Without_stimulation')
                              & (xcr_table[ch] == True)]

    healthy_samples = list(xcr_table.Path.apply(lambda x: x + file))
    # calculating fraction groups for healthy cohort

    h_f_extra_rare, h_f_rare, h_f_extra_small, h_f_small, h_f_med, h_f_large, h_f_extra_large, h_f_hyp, h_f_extra_hyp = calc_fraction_groups(
        healthy_samples, chain)

    h_f_extra_rare = statistics.mean(h_f_extra_rare)
    h_f_rare = statistics.mean(h_f_rare)
    h_f_small = statistics.mean(h_f_small)
    h_f_extra_small = statistics.mean(h_f_extra_small)
    h_f_med = statistics.mean(h_f_med)
    h_f_large = statistics.mean(h_f_large)
    h_f_extra_large = statistics.mean(h_f_extra_large)
    h_f_hyp = statistics.mean(h_f_hyp)
    h_f_extra_hyp = statistics.mean(h_f_extra_hyp)
    return h_f_extra_rare, h_f_rare, h_f_extra_small, h_f_small, h_f_med, h_f_large, h_f_extra_large, h_f_hyp, h_f_extra_hyp

def fraction_barplot(samples,
                     names,
                     chain,
                     material_type,
                     save_path,
                     healthy=False,
                     healthy_base='Target_Lab',
                     seq_type='Target',
                     xticklabels=True,
                     yticklabels=True,
                     x_anchor=1.6,
                     DPI=150,
                     legend=True,
                     fontsize=40,
                     figsize=(10, 10),
                     edgecolor='lightgrey'):
    """
    
    Plot fraction barplots for required samples.
    :param samples: list, file names
    :param names: list, sample names
    :param chain: str, required chain, can contain: alpha, beta, gamma, delta,
                  heavy, kappa, lambda, light, IgM, IgG
    :param material_type: str, tissue type: PBMC, T cells...
    :param healthy: bool, add "Healthy donor" column in barplot if True
    :param healthy_base: str, base for healthy cohort. 'Target' becouse function for target.
    :return : fig, ax
    
    """

    if type(samples) != list:
        samples = [samples]
        names = [names]
    # calculation fraction groups for all samples
    f_extra_rare, f_rare, f_extra_small, f_small, f_med, f_large, f_extra_large, f_hyp, f_extra_hyp = calc_fraction_groups(
        samples, chain)
    # if 'healthy' == True calculation fraction groups for healthy cohort
    if healthy:
        names.append('Healthy donors')
        h_f_extra_rare, h_f_rare, h_f_extra_small, h_f_small, h_f_med, h_f_large, h_f_extra_large, h_f_hyp, h_f_extra_hyp = fraction_healthy_cohort(
            chain, material_type, healthy_base, seq_type)

        f_extra_rare.append(h_f_extra_rare)
        f_rare.append(h_f_rare)
        f_extra_small.append(h_f_extra_small)
        f_small.append(h_f_small)
        f_med.append(h_f_med)
        f_large.append(h_f_large)
        f_extra_large.append(h_f_extra_large)
        f_hyp.append(h_f_hyp)
        f_extra_hyp.append(h_f_extra_hyp)

    # crate fraction groups dataframe
    fr_res_df = pd.DataFrame(
        {
            '10⁻¹ < Extra Hyperexpanded': f_extra_hyp,
            '10⁻² < Hyperexpanded ≤ 10⁻¹': f_hyp,
            '10⁻³ < Extra Large ≤ 10⁻²': f_extra_large,
            '10⁻⁴ < Large ≤ 10⁻³': f_large,
            '10⁻⁵ < Medium ≤ 10⁻⁴': f_med,
            '10⁻⁶ < Small ≤ 10⁻⁵': f_small,
            '10⁻⁷ < Extra Small ≤ 10⁻⁶': f_extra_small,
            '10⁻⁸ < Rare ≤ 10⁻⁷': f_rare,
            'Extra Rare ≤ 10⁻⁸': f_extra_rare
        },
        index=names)
    # set colour pallet
    hl_palette = {
        'Extra Rare ≤ 10⁻⁸': '#B40426',
        '10⁻⁸ < Rare ≤ 10⁻⁷': '#FF4E3C',
        '10⁻⁷ < Extra Small ≤ 10⁻⁶': '#FF8166',
        '10⁻⁶ < Small ≤ 10⁻⁵': '#F5B946',
        '10⁻⁵ < Medium ≤ 10⁻⁴': '#17C98D',
        '10⁻⁴ < Large ≤ 10⁻³': '#AFD7EF',
        '10⁻³ < Extra Large ≤ 10⁻²': '#88AADD',
        '10⁻² < Hyperexpanded ≤ 10⁻¹': '#6E7096',
        '10⁻¹ < Extra Hyperexpanded': '#394A59'
    }

    # image plotting
    fig, ax = plt.subplots(1, 1, figsize=figsize, facecolor='white', dpi=DPI)
    fr_res_df.plot.bar(ax=ax,
                       stacked=True,
                       color=hl_palette,
                       grid=False,
                       edgecolor='white',
                       width=0.9,
                       legend=False)
    if legend:
        ax.legend(loc='upper center',
                  bbox_to_anchor=(x_anchor, 1.0),
                  ncol=1,
                  fontsize=fontsize,
                  edgecolor=edgecolor)
    if 'Ig' in chain:
        ax.set_title(material_type + ' ' + chain, fontsize=fontsize)
    else:
        ax.set_title(material_type + ' ' + chain.capitalize(),
                     fontsize=fontsize)
    if xticklabels:
        ax.set_xticklabels(ax.get_xticklabels(),
                           rotation=60,
                           ha='right',
                           fontsize=fontsize)
    else:
        ax.set_xticklabels([], fontsize=fontsize)
        ax.get_xaxis().set_ticks([])
    if yticklabels:
        ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2],
                           fontsize=fontsize)
    else:
        ax.set_yticklabels([], fontsize=fontsize)
        ax.get_yaxis().set_ticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor(edgecolor)
#     return fig, ax
    plt.savefig('{}fraction_barplot_{}.svg'.format(save_path, chain), bbox_inches='tight')

def single_sample_fraction_barplots(samples,
                                    names,
                                    save_path,
                                    xcr,
                                    seq_type,
                                    material_type,
                                    base='Lab'):
    if xcr == 'bcr':
        chains = ['heavy', 'kappa', 'lambda', 'light']
        xtick_list = [True, True, True, True]
        ytick_list = [True, True, False, False]
        legend_list = [False, False, True, True]
    elif xcr == 'tcr-ab':
        chains = ['alpha', 'beta']
        xtick_list = [True, True]
        ytick_list = [True, False]
        legend_list = [False, True]
    else:
        chains = ['gamma', 'delta']
        xtick_list = [True, False]
        ytick_list = [True, True]
        legend_list = [False, True]
    for i in range(len(chains)):
        chain = chains[i]
        xticklebels = xtick_list[i]
        yticklebels = ytick_list[i]
        legend = legend_list[i]
        fraction_barplot(samples,
                         names,
                         chain,
                         material_type,
                         save_path,
                         healthy=True,
                         healthy_base='Lab',
                         seq_type=seq_type,
                         xticklabels=xticklebels,
                         yticklabels=yticklebels,
                         x_anchor=1.6,
                         DPI=150,
                         legend=legend,
                         fontsize=40,
                         figsize=(10, 10),
                         edgecolor='lightgrey')

# if __name__ == "__main__":
#     # creating 4 barplots for 2x2 table for 2 chains and 2 material_types
#     material_types = ['PBMC', 'B cells']
#     ch = ['heavy', 'light']
#     samples_1 = [
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_32/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_33/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_36/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_37/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_40/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_41/BCR_qc_pass.txt'
#     ]
#     samples_2 = [
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_30/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_31/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_34/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_35/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_38/BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_39/BCR_qc_pass.txt'
#     ]
#     l_names = [['18/18', '18/18', '18/16', '18/16', '18/15', '18/15'],
#                ['18/18', '18/18', '18/16', '18/16', '18/15', '18/15']]
#     healthy = True

#     for i in range(len(material_types)):
#         material_type = material_types[i]
#         names = l_names[i]
#         if i == 0:
#             samples = samples_1
#         else:
#             samples = samples_2
#         for j in range(len(ch)):
#             chain = ch[j]
#             if i == 0 and j == 0:
#                 yticklabels = True
#                 xticklabels = False
#                 legend = False
#             elif i == 0 and j == 1:
#                 yticklabels = True
#                 xticklabels = True
#                 legend = False
#             elif i == 1 and j == 0:
#                 yticklabels = False
#                 xticklabels = False
#                 legend = True
#             elif i == 1 and j == 1:
#                 yticklabels = False
#                 xticklabels = True
#                 legend = False
#             fraction_barplot(samples,
#                              names,
#                              chain,
#                              material_type,
#                              healthy=healthy,
#                              healthy_base='Target_Lab',
#                              seq_type='Target',
#                              xticklabels=xticklabels,
#                              yticklabels=yticklabels,
#                              x_anchor=1.6,
#                              DPI=150,
#                              legend=legend,
#                              fontsize=40,
#                              figsize=(10, 10),
#                              edgecolor='lightgrey')
#             if healthy == True:
#                 names.remove('Healthy donors')