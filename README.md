# Volcano Plot in Python

Inspired by R package [kevinblighe/EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano)

![image](https://user-images.githubusercontent.com/70199914/182571091-b04b7881-e5a1-4797-ba82-e0e2f1aedab7.png)


Input: Pandas df

  ![image](https://user-images.githubusercontent.com/70199914/182568630-4c9ab4cf-8afc-46f4-b51c-1580005803e3.png)



```
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from adjustText import adjust_text

def volcanoplot(res=rank, pval_cutoff=0.95, pval_colour_threshold=0.05, log2FC_colour_threshold=1, pval_label_cutoff=0.05, log2FC_label_cutoff=0.5, dotsize=4, title='Volcano Plot'):

    res=rank # Pandas df
    pval_cutoff=pval_cutoff # exclude all genes > pval_cutoff, as these swamp the plot
    pval_colour_threshold=pval_colour_threshold # threshold for colouring dots
    log2FC_colour_threshold=log2FC_colour_threshold # threshold for colouring dots
    pval_label_cutoff=pval_label_cutoff  # cutoff for dot labels
    log2FC_label_cutoff=log2FC_label_cutoff  # cutoff for dot labels
    dotsize=dotsize
    title=title
    
    toplot = res[res.pvals_adj <= pval_cutoff]

    # plot
    # plot non-significant genes with small log2FC
    plt.plot(toplot[(toplot.pvals_adj > pval_colour_threshold) & (toplot['log2FC'].abs()<log2FC_colour_threshold)].log2FC,
             toplot[(toplot.pvals_adj > pval_colour_threshold) & (toplot['log2FC'].abs()<log2FC_colour_threshold)].nlog10_pval_adj, 'o',
             color='#808080', alpha=.6, ms=dotsize, label='NS') # green

    # plot non-significant genes with large log2FC
    plt.plot(toplot[(toplot.pvals_adj > pval_colour_threshold) & (toplot['log2FC'].abs()>=log2FC_colour_threshold)].log2FC,
             toplot[(toplot.pvals_adj > pval_colour_threshold) & (toplot['log2FC'].abs()>=log2FC_colour_threshold)].nlog10_pval_adj, 'o',
             color='#1a9641', alpha=.6, ms=dotsize, label='log2FC') # grey 

    # plot significant genes with small log2FC
    plt.plot(toplot[(toplot.pvals_adj<=pval_colour_threshold) & (toplot['log2FC'].abs()<log2FC_colour_threshold)].log2FC,
             toplot[(toplot.pvals_adj<=pval_colour_threshold) & (toplot['log2FC'].abs()<log2FC_colour_threshold)].nlog10_pval_adj, 'o',
             color='#6495ED', alpha=.6, ms=dotsize, label='p-value') # blue

    # plot significant genes with large log2FC
    plt.plot(toplot[(toplot.pvals_adj<=pval_colour_threshold) & (toplot['log2FC'].abs()>=log2FC_colour_threshold)].log2FC,
             toplot[(toplot.pvals_adj<=pval_colour_threshold) & (toplot['log2FC'].abs()>=log2FC_colour_threshold)].nlog10_pval_adj, 'o',
             color='#FF3131', alpha=.6, ms=dotsize, label='p-value & log2FC') # red

    # axis labels etc
    plt.xlabel('log2FC')
    plt.ylabel('-log10(p)')
    plt.title(title)
    plt.legend(frameon=True, fontsize=12)

    # dot labels
    main_x = toplot[(toplot.pvals_adj<=pval_label_cutoff) & (toplot['log2FC'].abs()>=log2FC_label_cutoff)].log2FC
    main_y = toplot[(toplot.pvals_adj<=pval_label_cutoff) & (toplot['log2FC'].abs()>=log2FC_label_cutoff)].nlog10_pval_adj

    texts = []
    for x, y, s in zip(main_x, main_y, list(main_x.index)):
        texts.append(plt.text(x, y, s))

    adjust_text(texts,force_text=(1,1),arrowprops=dict(arrowstyle="-",lw=1))
    
    return(plt)

```

