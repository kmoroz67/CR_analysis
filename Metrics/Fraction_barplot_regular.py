import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics
import warnings
warnings.simplefilter('ignore')

PATH_TO_XCR_TABLE = '/uftp/Blood/prototype/local_databases/XCR_exist_chain_04_03_22.txt'

def calc_fraction_groups(samples, chain):
    """
    Calc Hyperexpanded, Large, Medium, Small, Rare fraction groups for required samples.
    :param samples: list, file names
    :param chain: str, required chain, can contain: alpha, beta, gamma, delta,
                  heavy, kappa, lambda, light, IgM, IgG
    :return : f_rare, f_small, f_med, f_large, f_hyp
    """

    fractions_list = []
    # data collection for each sample
    for samp in samples:
        freqs = {
            'Rare': 0,
            'Small': 0,
            'Medium': 0,
            'Large': 0,
            'Hyperexpanded': 0
        }

        # loading table and chain restriction
        if type(samp) == str:
            data = pd.read_csv(samp, sep='\t')
        else:
            data = samp

        if chain == 'light':
            data = data[(data.chain == 'lambda')
                        | (data.chain == 'kappa')].reset_index(drop=True)
        elif 'Ig' in chain:
            data = data[data.isotype == chain].reset_index(drop=True)
        else:
            data = data[data.chain == chain].reset_index(drop=True)

        # fraction calculation if chain is 'light', 'IgG', 'IgM'
        if chain == 'light' or 'Ig' in chain:
            s = data.clonotype_count.sum()
            data.fraction = data.clonotype_count.apply(lambda x: x / s)

        # data separating by fraction
        freqs['Rare'] += data[data.fraction <= 0.00001].fraction.sum()
        freqs['Small'] += data[(data.fraction <= 0.0001)
                            & (data.fraction > 0.00001)].fraction.sum()
        freqs['Medium'] += data[(data.fraction <= 0.001)
                             & (data.fraction > 0.0001)].fraction.sum()
        freqs['Large'] += data[(data.fraction <= 0.01)
                            & (data.fraction > 0.001)].fraction.sum()
        freqs['Hyperexpanded'] += data[data.fraction > 0.01].fraction.sum()

        fractions_list.append(freqs)

    f_rare = [samp['Rare'] for samp in fractions_list]
    f_small = [samp['Small'] for samp in fractions_list]
    f_med = [samp['Medium'] for samp in fractions_list]
    f_large = [samp['Large'] for samp in fractions_list]
    f_hyp = [samp['Hyperexpanded'] for samp in fractions_list]
    return f_rare, f_small, f_med, f_large, f_hyp

def fraction_healthy_cohort(chain, material_type, healthy_base, seq_type, age):
    '''
    Calc Hyperexpanded, Large, Medium, Small, Rare groups for healthy cohort (Target).
    :param chain: str, required chain, can contain: alpha, beta, gamma, delta,
                  heavy, kappa, lambda, light, IgM, IgG
    :param material_type: str, tissue type: PBMC or whatever you want (for 'Target_Lab'
                          if not PBMC, then T_cells or B_cells samples will be used).
    :param healthy_base: str, base for healthy cohort.
    :return : h_f_rare, h_f_small, h_f_med, h_f_large, h_f_hyp
    '''

    # loading path table
    xcr_table = pd.read_csv(PATH_TO_XCR_TABLE, sep='\t')
    if age !=None:
        xcr_table = xcr_table[xcr_table['Age'].apply(lambda x: True in [
        y in x for y in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
        ])].reset_index(drop=True)
        xcr_table['Age'] = xcr_table['Age'].astype(float)
        xcr_table = xcr_table.query('Age>({}-10)&Age<({}+10)'.format(age, age))
    if chain in ['alpha', 'beta', 'gamma', 'delta']:
        file = 'TCR_qc_pass.txt'
        ch = chain
        if healthy_base == 'Lab':
            if material_type == 'PBMC':
                t_m = ['PBMC']
            elif material_type == 'Whole_blood':
                t_m = ['PBMC', 'WBC']
            else:
                t_m = ['B_cells']
        else:
            t_m = [material_type]
    elif chain in ['heavy', 'kappa', 'lambda', 'light'] or 'Ig' in chain:
        file = 'BCR_qc_pass.txt'
        if 'Ig' in chain:
            ch = 'heavy'
        else:
            ch = chain
        if healthy_base == 'Lab':
            if material_type == 'PBMC':
                t_m = ['PBMC']
            elif material_type == 'Whole_blood':
                t_m = ['PBMC', 'WBC']
            else:
                t_m = ['B_cells']
        else:
            t_m = [material_type]
    else:
        print('ERROR: WRONG CHAIN')

    # path table restriction by 'healthy_base', 'chain', 'material_type'
    if chain=='light':
        xcr_table = xcr_table[(xcr_table.Base == healthy_base)
                              & (xcr_table.Cell_type.isin([x.replace(' ', '_') for x in t_m]))
                              & (xcr_table.Type == seq_type)
                              & (xcr_table.Diagnosis == 'Healthy')
                              &
                              (xcr_table.Stimulation == 'Without_stimulation')]
        xcr_table = xcr_table[xcr_table['lambda'] | xcr_table.kappa]
    else:
        xcr_table = xcr_table[(xcr_table.Base == healthy_base)
                              & (xcr_table.Cell_type.isin([x.replace(' ', '_') for x in t_m]))
                              & (xcr_table.Type == seq_type)
                              & (xcr_table.Diagnosis == 'Healthy')
                              &
                              (xcr_table.Stimulation == 'Without_stimulation')
                              & (xcr_table[ch] == True)]
    healthy_samples = list(xcr_table.Path.apply(lambda x: x + file))
    print(len(healthy_samples))
    # calculating fraction groups for healthy cohort
    h_f_rare, h_f_small, h_f_med, h_f_large, h_f_hyp = calc_fraction_groups(
        healthy_samples, chain)
    h_f_rare = statistics.mean(h_f_rare)
    h_f_small = statistics.mean(h_f_small)
    h_f_med = statistics.mean(h_f_med)
    h_f_large = statistics.mean(h_f_large)
    h_f_hyp = statistics.mean(h_f_hyp)
    return h_f_rare, h_f_small, h_f_med, h_f_large, h_f_hyp

def fraction_barplot(samples,
                     names,
                     chain,
                     material_type,
                     save_path,
                     age=None,
                     healthy=False,
                     healthy_base='Lab',
                     seq_type='BULK',
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
    f_rare, f_small, f_med, f_large, f_hyp = calc_fraction_groups(
        samples, chain)
    # if 'healthy' == True calculation fraction groups for healthy cohort
    if healthy:
        names.append('Healthy donors')
        h_f_rare, h_f_small, h_f_med, h_f_large, h_f_hyp = fraction_healthy_cohort(
            chain, material_type, healthy_base, seq_type, age)

        f_rare.append(h_f_rare)
        f_small.append(h_f_small)
        f_med.append(h_f_med)
        f_large.append(h_f_large)
        f_hyp.append(h_f_hyp)

    # crate fraction groups dataframe
    fr_res_df = pd.DataFrame(
        {
            '10⁻² < Hyperexpanded': f_hyp,
            '10⁻³ < Large ≤ 10⁻²': f_large,
            '10⁻⁴ < Medium ≤ 10⁻³': f_med,
            '10⁻⁵ < Small ≤ 10⁻⁴': f_small,
            'Rare ≤ 10⁻⁵': f_rare
        },
        index=names)
    # set colour pallet
    hl_palette = {
        'Rare ≤ 10⁻⁵': '#FF8166',
        '10⁻⁵ < Small ≤ 10⁻⁴': '#F4AE63',
        '10⁻⁴ < Medium ≤ 10⁻³': '#17C98D',
        '10⁻³ < Large ≤ 10⁻²': '#A6C4EA',
        '10⁻² < Hyperexpanded': '#6E7096'
    }

    # image plotting
    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=DPI)
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
                  fontsize=fontsize,
                  edgecolor=edgecolor)
    if 'Ig' in chain:
        ax.set_title(material_type + ' ' + chain, fontsize=fontsize)
    else:
        ax.set_title(material_type.replace('_',' ') + ' ' + chain.capitalize(),
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
    plt.tick_params(bottom=False, left=False)
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

if __name__ == "__main__":
    # creating 4 barplots for 2x2 table for 2 chains and 2 material_types
    material_types = ['PBMC', 'B cells']
    ch = ['heavy', 'light']
    samples_1 = [
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_32/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_33/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_36/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_37/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_40/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_41/BCR_qc_pass.txt'
    ]
    samples_2 = [
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_30/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_31/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_34/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_35/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_38/BCR_qc_pass.txt',
        '/uftp/Blood/db_calc_pipeline/target_lab/210708_NovaA_Sample_39/BCR_qc_pass.txt'
    ]
    l_names = [['18/18', '18/18', '18/16', '18/16', '18/15', '18/15'],
               ['18/18', '18/18', '18/16', '18/16', '18/15', '18/15']]
    healthy = True

    for i in range(len(material_types)):
        material_type = material_types[i]
        names = l_names[i]
        if i == 0:
            samples = samples_1
        else:
            samples = samples_2
        for j in range(len(ch)):
            chain = ch[j]
            if i == 0 and j == 0:
                yticklabels = True
                xticklabels = False
                legend = False
            elif i == 0 and j == 1:
                yticklabels = True
                xticklabels = True
                legend = False
            elif i == 1 and j == 0:
                yticklabels = False
                xticklabels = False
                legend = True
            elif i == 1 and j == 1:
                yticklabels = False
                xticklabels = True
                legend = False
            fraction_barplot(samples,
                             names,
                             chain,
                             material_type,
                             healthy=healthy,
                             healthy_base='Lab',
                             seq_type='BULK',
                             xticklabels=xticklabels,
                             yticklabels=yticklabels,
                             x_anchor=1.6,
                             DPI=150,
                             legend=legend,
                             fontsize=40,
                             figsize=(10, 10),
                             edgecolor='lightgrey')
            if healthy == True:
                names.remove('Healthy donors')