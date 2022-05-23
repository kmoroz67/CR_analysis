import matplotlib.pyplot as plt
import seaborn as sns

def top_clones_histogram(df,
                         top_n,
                         chain,
                         save_path,
                         material_type,
                         seq_type,
                         max_y_lim=None,
                         figsize=(10, 8),
                         fontsize=20,
                         DPI=150):

    if chain == 'light':
        df = df[df.chain.isin(['kappa', 'lambda'])]
        sum_i = df['clonotype_count'].sum()
        df['fraction'] = df['clonotype_count'].apply(lambda x:
                                                     (x / sum_i) * 100)
    elif 'Ig' in chain:
        df = df[df.isotype == chain]
        sum_i = df['clonotype_count'].sum()
        df['fraction'] = df['clonotype_count'].apply(lambda x:
                                                     (x / sum_i) * 100)
    else:
        df = df[df.chain == chain]
        df['fraction'] = df['fraction'].apply(lambda x: x * 100)

    df.sort_values(by='fraction', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df = df.head(top_n)
    if 'Ig' in chain:
        ch = chain
    else:
        ch = chain.capitalize()
    df['type'] = df['fraction'].apply(
        lambda x: '{} chain fraction ≥ 1%'.format(ch)
        if x >= 1 else '{} chain fraction < 1%'.format(ch))

    if chain in ['alpha', 'gamma', 'kappa', 'lambda', 'light']:
        color = 'tomato'
    else:
        color = 'steelblue'
    palette = {
        '{} chain fraction ≥ 1%'.format(ch): color,
        '{} chain fraction < 1%'.format(ch): '#E4E6EA'
    }
    fig, ax = plt.subplots(1, 1, facecolor='white', figsize=figsize, dpi=DPI)
    sns.barplot(data=df,
                x='cdr3aa',
                y='fraction',
                hue='type',
                palette=palette,
                dodge=False)
    if max_y_lim != None:
        ax.set_ylim([0, max_y_lim])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.grid(False)
    #     plt.legend([], [], frameon=False)
    plt.legend(prop={'size': fontsize})

    #     return fig, ax

    plt.savefig('{}{}_{}_top_{}_{}.svg'.format(save_path, material_type,
                                               seq_type, top_n, chain),
                bbox_inches='tight')