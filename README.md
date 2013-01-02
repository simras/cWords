

Example run:

jruby -J-Xmx4g scripts/cwords.rb -s <fasta-file> -r <rank-file> -w 6,7,8 -b 2

Options:

Usage: cwords [options]
    -w, --wordsize ARG               word length(s) you wish to search in (default 6,7)
    -b, --bg ARG                     Order of Markov background nucleotide model (default 2)
    -t, --threads ARG                use multiple threads to parallelize computations (default 1)
        --onlyanno                   only process annotated (i.e. mirbase) words
    -C, --custom_IDs                 Use your own sequences with matching IDs in rank and sequence files
    -A, --anders_ids                 Use Anders' integer IDs
    -x, --rank_split_mean            Split ranked list at mean
    -m, --rank_split_median          split ranked list at median
    -i, --rank_inverse               inverse all ranked lists
    -a, --rank_abs                   rank by absolute value
        --lead ARG                   Consider n most up regulated genes when looking for positively correlated words and visa versa for down regulated genes
        --inv_lead ARG               Consider most down regulated genes when looking for positively correlated words and visa versa for up regulated genes
        --interval ARG               Consider (like --interval n1,n2) an interval of genes (in the order of decreasing foldchange)
    -d, --rank_dist                  Analyse only leading genes with Z-score > 0.5 or Z-score < -0.5 in the respective up and down regulation analysis
    -r, --rankfile ARG               Rank file with IDs and one or more columns (will be mean collapsed) of a metric of expression changes (tab or space delimiter) or just one column with ordered IDs most down-regulated to most up-regulated
    -s, --seqfile ARG                Sequence file - Rank and sequence IDs should be one of the compatible IDs for the species you use
        --gen_plot ARG               Generate plot data and plots for the top k words, takes 2 passes.
    -u, --dump ARG                   Dump top words in file - 0 means all
        --report_words ARG           Report only on following words (comma separated)
        --plot_words ARG             Only make plots for following words (comma separated)
        --one_plot_words ARG         Make one plot with all following words (requires -x option to be set,comma separated)
        --test                       Testing mode - log files generated
    -R, --report_top ARG             Report top n words
    -g, --gene_set ARG               Print out what genes indicated words are enriched in. write words as comma-separated list.
        --plot_neg_words ARG         only make enrichment Enrichment Profile Plot files for negatively correlated words (comma separated)
        --plot_pos_words ARG         only make enrichment Enrichment Profile Plot files for positively correlated words (comma separated)
        --mkplot ARG                 Make Word Cluster Plot - highlight (with black border) the 8mer seed site (ex for miR-1: ACATTCCA) and its corresponding 7mer and 6mer seed sites. To highlight nothing write a word not in the [acgt] alphabet.
        --web                        Web server print mode.
    -N, --noAnn                      No miRNA-annotation on the Word Cluster Plot
        --annotFile ARG              Supply you own annotations, for Word Cluster Plot and word ranking.
        --species ARG                Different ID systems are used for different species, what's the species of your data? Currently we support Human Ensembl sequences as default (write human), Mouse (write mouse), Fruit Fly (write fruitfly) and Roundtype  Worm (write roundworm)

Results explained:

invalid IDs

positive vs negative set

header:
rank            word            z-score         p-value         fdr             ledge           annotation  