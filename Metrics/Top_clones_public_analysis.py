import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def top_clonotypes_public_analysis(df,
                                   chain,
                                   save_path,
                                   cohort_base,
                                   cohort_interest,
                                   top,
                                   material_type,
                                   seq_type,
                                   figsize=(10, 10),
                                   DPI=150):
    if chain in ['alpha', 'beta', 'gamma', 'delta']:
        xcr = 'TCR'
    else:
        xcr = 'BCR'
    if chain == 'heavy' or 'Ig' in chain:
        columns = ['cdr3aa', 'v_region', 'j_region', 'isotype']
    else:
        columns = ['cdr3aa', 'v_region', 'j_region']
    if 'Ig' in chain:
        ch = 'isotype'
    else:
        ch = 'chain'


#     if xcr == 'TCR':
#         df = df_tcr[df_tcr[ch] == chain].head(top)
#     else:
#         df = df_bcr[df_bcr[ch] == chain].head(top)
    df = df[df[ch] == chain].sort_values(by=['fraction'],
                                         ascending=False).head(top)
    if 'Ig' in chain:
        chain = 'heavy'
    df.set_index(columns, inplace=True, drop=True)
    df = df[['fraction']]
    columns = pd.MultiIndex.from_tuples([('Metrics', 'Fraction')],
                                        names=['Cohorts', 'Data'])
    df.set_axis(columns, axis='columns', inplace=True)
    public_df = pd.read_pickle(
        '/uftp/Blood/db_calc_pipeline/public_analysis/Test_cohorts_{}_aa_third_{}.pkl'
        .format(xcr, chain))
    public_df = public_df.loc[public_df.index.isin(df.index)]
    #  костыль
    if cohort_base == 'Regular':
        for col in ['Sample_counts', 'Cumulative_counts', 'Diagnoses_counts']:
            public_df[('Regular', col)] = public_df[[('Healthy', col),
                                                     ('Regular', col)
                                                     ]].apply(lambda x: sum(x),
                                                              axis=1)
    #
    df_concat = pd.concat([public_df, df], axis=1, sort=False)
    df_concat.sort_values(by=[('Metrics', 'Fraction')],
                          ascending=False,
                          inplace=True)
    df_concat[('Metrics',
               '{}/{}'.format(cohort_interest, cohort_base))] = df_concat[[
                   (cohort_base, 'Sample_counts'),
                   (cohort_base, 'Cumulative_counts'),
                   (cohort_interest, 'Sample_counts'),
                   (cohort_interest, 'Cumulative_counts')
               ]].apply(lambda x: ((x[3] / x[2]) / (x[1] / x[0]))
                        if x[2] != 0 else 0,
                        axis=1)
    df_concat.fillna(0, inplace=True)
    top_clonotypes_public_analysis_plot(df_concat, save_path, chain,
                                        material_type, seq_type, cohort_base,
                                        cohort_interest, top, figsize, DPI)
                                        
def top_clonotypes_public_analysis_plot(df, save_path, chain, material_type,
                                        seq_type, cohort_base, cohort_interest,
                                        top, figsize, DPI):
    df[('Metrics', 'cdr3aa')] = [x[0] for x in df.index]
    df[('Metrics', 'Type')] = df[(
        'Metrics',
        'Fraction')].apply(lambda x: '{} chain fraction ≥ 1%'.format(
            chain) if x >= 0.01 else '{} chain fraction < 1%'.format(chain))

    if chain in ['alpha', 'gamma', 'kappa', 'lambda', 'light']:
        color = 'tomato'
    else:
        color = 'steelblue'
    palette = {
        '{} chain fraction ≥ 1%'.format(chain): color,
        '{} chain fraction < 1%'.format(chain): '#E4E6EA'
    }
    fig, ax = plt.subplots(2, 2, facecolor='white', figsize=figsize, dpi=DPI)
    sns.barplot(data=df,
                x=('Metrics', 'cdr3aa'),
                y=('Metrics', 'Fraction'),
                hue=('Metrics', 'Type'),
                ax=ax[0][0],
                palette=palette,
                dodge=False)
    sns.barplot(data=df,
                x=('Metrics', 'cdr3aa'),
                y=(cohort_base, 'Sample_counts'),
                ax=ax[1][0],
                dodge=False,
                log=True)
    sns.barplot(data=df,
                x=('Metrics', 'cdr3aa'),
                y=('Metrics', '{}/{}'.format(cohort_interest, cohort_base)),
                ax=ax[0][1],
                dodge=False)
    sns.barplot(data=df,
                x=('Metrics', 'cdr3aa'),
                y=(cohort_interest, 'Sample_counts'),
                ax=ax[1][1],
                dodge=False,
                log=True)
    ax[0][0].set_xlabel('')
    ax[0][0].set_ylabel('Fraction')
    ax[0][1].set_xlabel('')
    ax[0][1].set_ylabel('{}/{}'.format(cohort_interest, cohort_base))
    ax[0][0].set_xticklabels([])
    ax[0][1].set_xticklabels([])
    ax[0][0].set_xticks([])
    ax[0][1].set_xticks([])
    ax[1][0].set_xlabel('')
    ax[1][0].set_ylabel(cohort_base)
    ax[1][0].set_xticklabels(ax[1][0].get_xticklabels(), rotation=90)
    ax[1][1].set_xlabel('')
    ax[1][1].set_ylabel(cohort_interest)
    ax[1][1].set_xticklabels(ax[1][1].get_xticklabels(), rotation=90)
    
    handles, labels = ax[0][0].get_legend_handles_labels()
    ax[0][0].legend(handles=handles[1:], labels=labels[1:])
    
    
    fig.suptitle('{} {} Top {} clonotypes public analysis'.format(
        material_type, seq_type, top),
                 fontsize=18)
    plt.savefig('{}Public_{}_{}_top_{}_{}.svg'.format(save_path, material_type,
                                                      seq_type, top, chain),
                bbox_inches='tight')