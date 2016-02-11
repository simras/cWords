cWords is a software tool that measure correlations of short oligonucleotide sequence (word) occurrences and i.e. expression changes in a two condition experiment. In summary, it produces a statistic that quantifies overrepresented a word is near the extremity of a ranked list of sequences.

## 1. REQUIREMENTS #

All software components require Ruby (>=1.8.6, http://www.ruby-lang.org/) and JRuby (>=1.4.0, http://jruby.org/).
The software has only been tested on a Unix platform (both Linux and Mac OS X Lion).


## 2. INSTALL #

* Install Ruby (www.ruby-lang.org, check if it is already installed: 'ruby -v')
* Install Java (www.java.com, check if it is already installed: 'java -version'),
* Install JRruby (www.jruby.org, check if it is already installed: 'jruby -v'),
  make sure that you have the 'jruby' command in your path.

## 3. USAGE #
Below is a short summary of how to use the software. The full set of options for each script can be be listed from the command line by using the '-h' flag. The most important ones are the following:

### Options:

    Usage: cwords [options]
        -w, --wordsize ARG               word length(s) you wish to search in (default 6,7)
        -b, --bg ARG                     Order of Markov background nucleotide model (default 0)
        -t, --threads ARG                use multiple threads to parallelize computations (default 1)
        -C, --custom_IDs                 Use your own sequences with matching IDs in rank and sequence files
        -A, --anders_ids                 Use integer IDs
        -x, --rank_split_mean            Split ranked list at mean
        -r, --rankfile ARG               Rank file with IDs and one or more columns (will calculate mean across columns) of a metric of expression changes (tab or space delimiter) or just one column with ordered IDs most down-regulated to most up-regulated
        -s, --seqfile ARG                Sequence file - Rank and sequence IDs should be one of the compatible IDs for the species you use
            --gen_plot ARG               Generate plot data and plots for the top k words, takes 2 passes.
            --mkplot ARG                 Make Word Cluster Plot - highlight (with black border) the 8mer seed site (ex for miR-1: ACATTCCA) and its corresponding 7mer and 6mer seed sites. To highlight nothing write a word not in the [acgt] alphabet.
        -N, --noAnn                      No miRNA-annotation on the Word Cluster Plot
            --annotFile ARG              Supply you own annotations, for Word Cluster Plot and word ranking.
            --species ARG                Different ID systems are used for different species, what's the species of your data? Currently we support Human Ensembl sequences as default (write human), Mouse (write mouse), Fruit Fly (write fruitfly) and C. elegans (write roundworm)

 
### Example runs:
The following are standard analyses to run with cWords.
    
    tri-nucleotide background model, word sizes 6, 7 and 8, using 40 processors
    jruby -J-Xmx4g scripts/cwords_Mult_MEM.rb -s <fasta-file> -r <rank-file> -w 6,7,8 -b 2 -p 40

    IDs in you sequences and rank-file match
    jruby -J-Xmx4g scripts/cwords_Mult_MEM.rb -s <fasta-file> -r <rank-file> -w 6,7,8 -b 2 -C

    Mononucleotide background model, word cluster plot (highlighting ACATTCCA, ACATTCC, CATTCCA, ACATTC, CATTCC), Enrichment profile graphs for the 20 most significant words
    jruby -J-Xmx4g scripts/cwords_Mult_MEM.rb -s <fasta-file> -r <rank-file> -w 6,7,8 -b 0 -mkplot ACATTCCA -gen_plot 20

## 4. INTERPRETATION #
Results mainly compose of three elements. A ranking of most strongly correlated word, a word cluster plot and Enrichment profile plots. 

### Positive and negative set 
The analysis is divided into 2 passes; one where words that are over represented in up-regulated genes (ie. positively correlated words) and one where negatively correlated words is evalaluated. If all genes are considered in the analysis in each pass all words of the length in question will be divided into the negative set and positive set except the words that occur 5 or less times in the sequences. These sets can be divided in different ways (see options -h), and you can consider only most regulated genes in the two passes. If you consider the most down-regulated words in the negative pass amd the most up-regulated in the positive pass a word can be present in both the negative and the positive set.

### The output 
The final output of this analysis produces a summary of the top correlating words (a list for each end of the ranked list), i.e. words over-represented in beginning of list:

    Top 10 words
    rank            word            z-score         p-value         fdr             ledge           annotation     
    1               cactgcc         22.54           1.00e-10        6.58e-09        1651            hsa-miR-34a-5p,hsa-miR-34c-5p,hsa-miR-449a,hsa-miR-449b-5p
    2               actgcca         20.62           1.00e-10        6.98e-09        1115            hsa-miR-34a-5p,hsa-miR-34c-5p,hsa-miR-449a,hsa-miR-449b-5p,hsa-miR-548au-3p
    3               cctgccc         20.54           1.00e-10        6.03e-09        2702            hsa-miR-6721-5p
    4               ccctgcc         19.82           1.00e-10        6.17e-09        3252            hsa-miR-1207-5p,hsa-miR-4763-3p
    5               ccctggg         18.74           1.00e-10        6.22e-09        3500            hsa-miR-1915-3p
    6               ctgcccc         17.16           1.00e-10        5.73e-09        2546            hsa-miR-486-3p 
    7               ctgccct         16.97           1.00e-10        5.69e-09        3248            hsa-miR-4632-5p,hsa-miR-4436b-3p
    8               ccagccc         16.54           1.00e-10        6.47e-09        2690            hsa-miR-762,hsa-miR-4492,hsa-miR-4498,hsa-miR-5001-5p
    9               ccccagc         16.53           1.00e-10        6.32e-09        2706            hsa-miR-4731-5p
    10              ctgggcc         16.42           1.00e-10        5.90e-09        3524  
    ...

* 'z-score' is a correlation statistic for the given word after correction for correlations obtained from random gene list orderings.
* 'fdr' (false discovery rate) is the estimated proportion of false discoveries for the given z-score threshold.
* 'ledge' is the leading-edge which denotes the position in the gene list where the maximum imbalance was measured; genes before this threshold are relatively enriched for the word compared to genes after the threshold.

### Invalid IDs
If you do not use -C the IDs in the rank-file will be mapped to the sequences, when it was not possible to map the rank-file ID to a sequence we report the ID as invalid and. Problems with many invalid IDs typically occur when one uses ID associated with different versions or even assemblies (like sequences hg19 and microArray probes from hg18). Invalid IDs are written to a file facilitating further investigation (file name: invalid_ids.txt).

### Other concerns
Memory consumption and running time can vary a lot. The word length (-w) has a significant impact on this, if you want to run for word lengths > 8 works best with 4 GB memory or more (depending on number of genes in your analysis). Word lengths > 9 works best wih 10 GB or more.

Generally running time grows exponentially with word length and linearly with the number of bases in the sequences that need to be analysed.

## 4. How TO CITE
Rasmussen S, Jacobsen A and Krogh A;cWords - systematic microRNA regulatory motif discovery from mRNA expression data; Silence (2013)

## 5. LICENSE
Copyright (c) 2011, Simon H. Rasmussen.
The software is open source and released under the MIT license.
