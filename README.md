cWords is a software tool that measure correlations of short oligonucleotide sequence ('words') occurrences and i.e. expression changes in a given experiment. In summary, it produces a statistic that a given word is overrepresented near the extremity of a ranked list of sequences.

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
        -b, --bg ARG                     Order of Markov background nucleotide model (default 2)
        -t, --threads ARG                use multiple threads to parallelize computations (default 1)
        -C, --custom_IDs                 Use your own sequences with matching IDs in rank and sequence files
        -A, --anders_ids                 Use Anders' integer IDs
        -x, --rank_split_mean            Split ranked list at mean
        -r, --rankfile ARG               Rank file with IDs and one or more columns (will be mean collapsed) of a metric of expression changes (tab or space delimiter) or just one column with ordered IDs most down-regulated to most up-regulated
        -s, --seqfile ARG                Sequence file - Rank and sequence IDs should be one of the compatible IDs for the species you use
            --gen_plot ARG               Generate plot data and plots for the top k words, takes 2 passes.
            --mkplot ARG                 Make Word Cluster Plot - highlight (with black border) the 8mer seed site (ex for miR-1: ACATTCCA) and its corresponding 7mer and 6mer seed sites. To highlight nothing write a word not in the [acgt] alphabet.
        -N, --noAnn                      No miRNA-annotation on the Word Cluster Plot
            --annotFile ARG              Supply you own annotations, for Word Cluster Plot and word ranking.
            --species ARG                Different ID systems are used for different species, what's the species of your data? Currently we support Human Ensembl sequences as default (write human), Mouse (write mouse) and Fruit Fly (write fruitfly)

 
### Example runs:
The following are standard analyses to run with cWords.    
    tri-nucleotide background model, word sizes 6, 7 and 8, using 40 processors
    jruby -J-Xmx4g scripts/cwords.rb -s <fasta-file> -r <rank-file> -w 6,7,8 -b 2 -p 40

    IDs in you sequences and rank-file match
    jruby -J-Xmx4g scripts/cwords.rb -s <fasta-file> -r <rank-file> -w 6,7,8 -b 2 -C

    Mononucleotide background model, Word Cluster Plot (highlighting ACATTCCA, ACATTCC, CATTCCA, ACATTC, CATTCC), Enrichment profile graphs for the 20 most significant words
    jruby -J-Xmx4g scripts/cwords.rb -s <fasta-file> -r <rank-file> -w 6,7,8 -b 0 -mkplot ACATTCCA -gen_plot 20

## 4. INTERPRETATION #
Results mainly compose of three elements. A ranking of most strongly correlated word, a Word Cluster Plot and Enrichment profile plots. 

### Positive and negative set 
The analysis is divided into 2 passes; one where words that are over represented in up-regulated genes (ie. positively correlated words) and one where negatively correlated words is evalaluated. If all genes are considered in the analysis in each pass all words of the length in question will be divided into the negative set and positive set except the words that occur 5 or less times in the sequences. These sets can be divided in different ways (see options -h), and you can consider only most regulated genes in the two passes. If you consider the most down-regulated words in the negative pass amd the most up-regulated in the positive pass a word can be present in both the negative and the positive set.

### The output 
The final output of this analysis produces a summary of the top correlating words (a list for each end of the ranked list), i.e. words over-represented in beginning of list:

    Top 10 words
    rank      word      RS        z-score   p-value   fdr       ledge
    1         ttatc     13.92     4.00      2.00e-03  7.94e-01  449
    2         atccc     12.79     3.78      3.99e-03  6.59e-01  417
    3         gtaatc    8.24      3.75      2.00e-03  5.01e-01  280
    4         gtgaaa    6.77      3.68      5.99e-03  4.56e-01  137
    5         cctat     10.26     3.61      3.99e-03  4.45e-01  384
    6         tttatc    11.93     3.59      3.99e-03  3.88e-01  449
    ...

* 'z-score' is a correlation statistic for the given word after correction for correlations obtained from random gene list orderings.
* 'fdr' (false discovery rate) is the estimated proportion of false discoveries for the given z-score threshold.
* 'ledge' is the leading-edge which denotes the position in the gene list where the maximum imbalance was measured; genes before this threshold are relatively enriched for the word compared to genes after the threshold.

### Plots
Write about the plot ref. the paper

### invalid IDs
If you do not use -C the IDs in the rank-file will be mapped to the sequences, when it was not possible to map the rank-file ID to a sequence we report the ID as invalid and. Problems with many invalid IDs typically occur when one uses ID associated with different versions or even assemblies (like sequences hg19 and microArray probes from hg18). Invalid IDs are written to a file facilitating further investigation (file name: invalid_ids.txt).

### Other concerns
Memory consumption and running time can vary significantly. 

The word length (-w) has a significant impact on the, if you want to run for word lengths > 8 works best with 4 GB memory or more (depending on number of genes in your analysis). Word lengths > 9 works best wih 10 GB or more.

Generally running time grows exponentially with word length and linearly with the number of sequences that need to be analysed.

## 4. How TO CITE

## 5. LICENSE
Copyright (c) 2011, Simon H. Rasmussen.
The software is open source and released under the MIT license (license is included).