#!/usr/bin/env jruby
#
# cwords - Regulatory Motif Discovery tool
# 
# Simon H. Rasmussen and Anders Jacobsen
# Center for Bioinformatics - University of Copenhagen
########################################################

srcdir = File.dirname(__FILE__)
basedir = File.expand_path("..", srcdir) + "/"
libdir = File.expand_path("..",  srcdir) + '/lib/'
resdir = File.expand_path("..",  srcdir) + '/resources/'

$LOAD_PATH << libdir

require libdir + 'wordRS-lib.rb'
require 'progressbar'
require 'optparse'
require 'rubygems'
require 'thread'
require 'java'
require 'BinomEvaluatorMult.jar'
require 'PermutationStats.jar'
require 'JavaFreeMem.jar'
java_import 'BinomEvaluatorMult'
java_import 'JavaFreeMem'
java_import 'PermutationStats'

# Default options
options = Hash.new
options[:wordsize] = [6,7]
options[:split_words] = nil
options[:scoring_scheme] = 'bpval'
options[:permutations] = 0
options[:seqshuffles] = 100
options[:rankfile] = nil
options[:seqfile] = nil
options[:report_words] = nil
options[:onlyanno] = nil
options[:dump] = 0
options[:testing] = nil
options[:rank_all] = true
options[:rank_inverse]=nil
options[:rank_split_median]=nil
options[:rank_abs]=nil
options[:gene_set]=nil
options[:lead]=nil
options[:inv_lead]=nil
options[:web]=false
options[:bg] = "1"
options[:threads]=1
options[:genplot]=0
options[:report_top]=10
options[:custom_IDs] = nil
options[:plot_neg_words]=[]
options[:plot_pos_words]=[]
options[:plot_words]=[]
options[:one_plot_words]=[]
options[:intv]=nil
options[:mkplot]=nil
options[:noAnn]=false
options[:annotFile] = ""
# NEW
options[:species] = "human"
options[:anders_ids] = false
options[:pa_seq] = false
options[:flatBgModel] = false

posTopWords = Array.new(options[:genplot])
negTopWords = Array.new(options[:genplot])

$coptions = OptionParser.new do |opts|
  opts.banner = "Usage: cwords [options]"
  
# analysis settings
  opts.on("-w", "--wordsize ARG", "word length(s) you wish to search in (default 6,7)") { |o| options[:wordsize] = o.split(",").map{|x| x.to_i}}
  opts.on("-b", "--bg ARG", "Order of Markov background nucleotide model (default 0)") {|o| options[:bg] = (o.to_i + 1).to_s}
  opts.on("-t", "--threads ARG", "use multiple threads to parallelize computations (default 1)") {|o| options[:threads] = o.to_i}
  opts.on(      "--onlyanno", "only process annotated (i.e. mirbase) words") {|o| options[:onlyanno] = true}
  opts.on("-C", "--custom_IDs", "Use your own sequences with matching IDs in rank and sequence files") {|o| options[:custom_IDs] = true}
  opts.on("-A", "--anders_ids", "Use Anders' integer IDs") {|o| options[:anders_ids] = true}

  # rank control
  opts.on("-x", "--rank_split_mean", "Split ranked list at mean") {|o| options[:rank_all] = false}
  opts.on("-m", "--rank_split_median", "split ranked list at median") {|o| options[:rank_split_median] = true}
  opts.on("-i", "--rank_inverse", "inverse all ranked lists") {|o| options[:rank_inverse] = true}
  opts.on("-a", "--rank_abs", "rank by absolute value") {|o| options[:rank_abs] = true}
  opts.on(      "--lead ARG", "Consider n most up regulated genes when looking for positively correlated words and visa versa for down regulated genes") {|o| options[:lead] = o.to_i}
  opts.on(      "--inv_lead ARG", "Consider most down regulated genes when looking for positively correlated words and visa versa for up regulated genes") {|o| options[:inv_lead] = o.to_i}
  opts.on(      "--interval ARG", "Consider (like --interval n1,n2) an interval of genes (in the order of decreasing foldchange)") {|o| options[:intv] = o.split(",")}
  opts.on("-d", "--rank_dist", "Analyse only leading genes with Z-score > 0.5 or Z-score < -0.5 in the respective up and down regulation analysis") {|o| options[:dist] = true}
  opts.on("-F", "--flat_bg", "Use a background model that assumes all nucleotides (and di or trinucleotides) are equally probable.") {|o| options[:flatBgModel] = true}

  # files and directories
  opts.on("-r", "--rankfile ARG", "Rank file with IDs and one or more columns (will be mean collapsed) of a metric of expression changes (tab or space delimiter) or just one column with ordered IDs most down-regulated to most up-regulated") {|o| options[:rankfile] = o}
  opts.on("-s", "--seqfile ARG", "Sequence file - Rank and sequence IDs should be one of the compatible IDs for the species you use") {|o| options[:seqfile] = o}
  opts.on("-P", "--PASeq", "By use of PolyA sequencing data you can assign the sequences a weight and the most probable UTR isoform will be selected.") {|o| options[:pa_seq] = true}

  # output control
  opts.on(      "--gen_plot ARG", "Generate plot data and plots for the top k words, takes 2 passes.") {|o| options[:genplot] = o.to_i }
  opts.on("-u", "--dump ARG", "Dump top words in file - 0 means all") { |o| options[:dump] = o.to_i}
  opts.on(      "--report_words ARG", "Report only on following words (comma separated)") {|o| options[:report_words] = o.split(',')}
  opts.on(      "--plot_words ARG", "Only make plots for following words (comma separated)") {|o| options[:plot_words] = o.split(',')}
  opts.on(      "--one_plot_words ARG", "Make one plot with all following words (requires -x option to be set,comma separated)") {|o| options[:one_plot_words] = o.split(',')}
  opts.on(      "--cooccur ARG", "Do following words cooccur") {|o| options[:cooccur] = o.split(',')}
  opts.on(      "--test", "Testing mode - log files generated") {|o| options[:testing] = true}
  opts.on("-R", "--report_top ARG", "Report top n words"){ |o| options[:report_top] = o.to_i}
  opts.on("-g", "--gene_set ARG", "Print out what genes indicated words are enriched in. write words as comma-separated list.") {|o| options[:gene_set] = o.split(',')}
  opts.on(      "--plot_neg_words ARG", "only make enrichment Enrichment Profile Plot files for negatively correlated words (comma separated)") {|o| options[:plot_neg_words] = o.split(',')}
  opts.on(      "--plot_pos_words ARG", "only make enrichment Enrichment Profile Plot files for positively correlated words (comma separated)") {|o| options[:plot_pos_words] = o.split(',')}
  opts.on(      "--mkplot ARG", "Make Word Cluster Plot - highlight (with black border) the 8mer seed site (ex for miR-1: ACATTCCA) and its corresponding 7mer and 6mer seed sites. To highlight nothing write a word not in the [acgt] alphabet.") {|o| options[:mkplot] = o.downcase}
  opts.on(      "--web", "Web server print mode.") {|o| options[:web] = true}
  opts.on("-N", "--noAnn", "No miRNA-annotation on the Word Cluster Plot") {|o| options[:noAnn] = true}
  opts.on("--annotFile ARG", "Supply you own annotations, for Word Cluster Plot and word ranking.") {|o| options[:annotFile] = o}
# NEW
  opts.on("-P", "--PASeq", "By use of PolyA sequencing data you can assign the sequences a weight and the most probable UTR isoform will be selected.") {|o| options[:pa_seq] = true}
  opts.on(      "--species ARG", "Different ID systems are used for different species, what's the species of your data? Currently we support Human Ensembl sequences as default (write human), Mouse (write mouse), Fruit Fly (write fruitfly) and Roundtype  Worm (write roundworm)") {|o| options[:species] = o.downcase}
end

# Parse inputs and show help if incomplete
def show_help(msg="", code=0, io=STDOUT)
  io.puts "#{msg}\n#{$coptions}"
  exit(code)
end

begin
  $coptions.parse!(ARGV) 
rescue OptionParser::ParseError => error
  puts error.message
  puts $coptions 
  exit
end

[:rankfile].each{|p| show_help("option '#{p}' mandatory") if options[p].nil?}
show_help("seqfile required") if !options[:seqfile]
show_help("scoring scheme must be: bpval") if !(['bpval'].include?(options[:scoring_scheme]))


# Gene intervals are converted for easier access
if options[:intv] then
  if options[:intv].length == 2 then
    gene_intv = Array(2)
    options[:intv].each_with_index{ |it,i|
      gene_intv[i] = it.to_i
    }
  end
end

# resolve conflict between options -A and -C
if options[:anders_ids] && options[:custom_IDs] then
  options[:custom_IDs] = false
end

# Checking if rank-file is correct format
# NEW
#print "egrep -x \'^[A-Za-z0-9.,_:\\+\;-]+\'","\n"
#print "egrep -x \'^[A-Za-z0-9.,_:\\+\;-]+[\\t ]+-?[0-9.]+\'","\n"
 
if system("head -1 " + options[:rankfile] + "|egrep -x \'^[A-Za-z0-9.,_:\\+\;-]+\' 1> /dev/null") then
  ok = true
  # put extra line
  system(libdir + "addLines.sh " + options[:rankfile])
  options[:rankfile] =  options[:rankfile] + ".tmp"
elsif system("head -1 " + options[:rankfile] + "|egrep -x \'^[A-Za-z0-9.,_:\\+\;-]+[\t ]+-?[0-9.]+\' 1> /dev/null") then
  ok = true
  # do nothing
else
  ok = false
end
raise ArgumentError, 'rank-file (-r) should have only one column with ordered IDs of tap-separate columns. ID can only be composed of following characters a-z,A-Z,0-9 and .,_-;:+' unless ok

testing = options[:testing]
output_top = options[:report_top]

# get filename without directory
rankfilename = File.basename(options[:rankfile])

# random string suffix of results folder
ranStr = random_string()
outDir = Dir.pwd() + "/" + rankfilename + "_" + ranStr + "/"

# species
if options[:species] == "mouse" then
  spepatt = "mmu"
elsif options[:species] == "roundworm" then
  spepatt = "cel"
elsif options[:species] == "fruitfly" then
  spepatt = "dme"
elsif options[:species] == "human" then
  spepatt = "hsa"
else
  raise "Sorry we don\'t have sequence ID infor for " + options[:species] + ", if we should please check your spelling, use only lower case letters or contact us."
end

if options[:annotFile] == "" then
  annofile = resdir  + "annotations_" + options[:species] + ".tsv"
else
  annofile = options[:annotFile]
end

# Disable miR-annotations
if options[:noAnn] && options[:annotFile] == "" then
  ann = "\"\""
else
  ann = annofile
end

# generate annotation file
if Dir[annofile][0].nil? then
  system("python " + resdir + "add_annotation.py -a " + resdir + "mature.fa -o " + resdir + "annotations_" + options[:species] + ".tsv -p " + spepatt)
end
#annofile = basedir + "/resources/" + "annot_gap2.tsv" #annotation

# read in annotations
word_annotation = Hash.new("") # seq => family
if !options[:noAnn] then
  IO.readlines(annofile).each{|l|word_annotation[l.split("\t")[0]] = l.split("\t")[1]}
end

# gene id mapping
if options[:species] == "human" then
  tidfile = basedir + "/resources/" + "genemap.tsv"
else
  tidfile = basedir + "/resources/" + "genemap_" + spepatt + ".csv"
end

seqshuffles = 5000

sequences = nil
if options[:bg] == "var" then
  varbg = true
else
  varbg = false
  bg = options[:bg].to_i
end

threads=options[:threads]

puts ">> Parameters"
options.each{|k,v| puts sprintf("%-20s:  %s",k,v) if !v.nil?} if !options[:web]


###
### Main program
###

# make output directory
system("mkdir " + outDir + " 1> /dev/null"+ " 2> /dev/null")

if testing then
  # make logfile for web-debugging pursposes
  logfile = File.new(outDir + rankfilename + "_log.txt","w")  
end

# web output
if options[:web]
  $stderr.puts outDir
end

# Make as default 1 pass. Second pass for plotting 
pass = true
newHash = Hash.new()
invFile = File.new(outDir + "invalid_ids.txt","w")

if options[:custom_IDs]
  # Custom_IDs map IDs from rank to sequence file, (seq/rank) data will only be included when there are matching IDs.
  seqmap = Hash.new
  idcounter = 0
  id_num = 0
  filtered = 0
  idmap = Hash.new
  allh = Hash.new {|h,k| h[k] = []}  
  IO.readlines(options[:rankfile]).each do |l|
    if /\t/ =~ l then
      l = l.split("\t")
    elsif / / =~ l then
      l = l.split(" ")
    end
    # map rank-IDs to integers
    idmap[l[0].chomp] = id_num.to_s
    newl=Array.new
    
    l[1..(l.size)].each do |val|
      newl << (val.chomp).to_f
    end
    allh[id_num.to_s] = newl
    id_num = id_num + 1
  end
  internal_ids = idmap.invert

  # Choose sequences that map to rank-file, if ID is ambigous choose longest sequence
  puts "\n>> reading sequences ..." if !options[:web] 
  sequences = Hash.new
#<<<<<<< HEAD
  num_reads = Hash.new if options[:pa_seq]
#=======
#  num_reads = Hash.new
#  pa_reads = false
#  nrc = 0
#>>>>>>> 6476e441b2297c9517df809724fb4e656904f90e
  IO.readlines(options[:seqfile],">").each do |entry|
    ls = entry.split("\n").map{|x| x.chomp}
    if ls[1..(ls.size)].join('') == "Sequence unavailable>" then
      next
    end
    # NEW
#<<<<<<< HEAD
    if options[:pa_seq] then
      tmpl = ls[0].split("|")
      ls[0] = tmpl[0]
      tmp_num_reads =  tmpl[1]
#=======
    #if options[:pa_seq] then
#    tmpl = ls[0].split("#")
#    print "--------- sharp -------------"
#    if tmpl.size > 1 then
#      pa_reads = true
#      ls[0] = tmpl[0]
#      tmp_num_reads =  tmpl[1].to_i
      #print tmp_num_reads, "\n"
#>>>>>>> 6476e441b2297c9517df809724fb4e656904f90e
    else
      tmpl = ls[0].split("|")
      if tmpl.size > 1 then
        # Biomart ENSEMBL format
        # 1st ID Genome
        #ls[0] = tmpl[0]
        # 2nd ID Transcript
        ls[0] = tmpl[1]
      end
    end
    if idmap[ls[0].chomp] then
      if sequences.has_key?(idmap[ls[0].chomp]) then
        newseq = ls[1..(ls.size)].join('').upcase.gsub('U','T').gsub('>','').gsub("SEQTENCE TNAVAILABLE",'')
#<<<<<<< HEAD
        # filter by length
        if options[:pa_seq]
          if tmp_num_reads > num_reads[idmap[ls[0].chomp]] then
#=======
       # filter by length
#        if pa_reads then
#          if tmp_num_reads > num_reads[idmap[ls[0].chomp]] then
            #print ls[0], " ",tmp_num_reads," ",newseq.size, "\n"
#            nrc = nrc + 1
#>>>>>>> 6476e441b2297c9517df809724fb4e656904f90e
            sequences[idmap[ls[0].chomp]] = newseq
          elsif tmp_num_reads == num_reads[idmap[ls[0].chomp]] && newseq.size > sequences[idmap[ls[0].chomp]].size then
            sequences[idmap[ls[0].chomp]] = newseq
          end
        else
          if newseq.size > sequences[idmap[ls[0].chomp]].size then
            sequences[idmap[ls[0].chomp]] = newseq
          end
        end
      else
        sequences[idmap[ls[0].chomp]] = ls[1..(ls.size)].join('').upcase.gsub('U','T').gsub('>','').gsub("SEQTENCE TNAVAILABLE",'') # last field is ">"
#<<<<<<< HEAD
        num_reads[idmap[ls[0].chomp]] = tmp_num_reads if options[:pa_seq]
#=======
#        num_reads[idmap[ls[0].chomp]] = tmp_num_reads if pa_reads
#>>>>>>> 6476e441b2297c9517df809724fb4e656904f90e
        if sequences[idmap[ls[0].chomp]].size <= 0 then
          sequences.delete(idmap[ls[0].chomp])
        end
      end
    end
  end
  #print "nrc ",nrc if !options[:web]
elsif options[:anders_ids] && options[:species] == "human" then
  # Use Anders' original ID framework - only Human
  filtered = 0  
  idmap = Hash.new
  internal_ids = Hash.new
  
  IO.readlines(tidfile).each do |l|
    ll = l.split(" ")
    tid = ll[0]
    ll[1].split(",").each{|extid| 
      idmap[extid] = tid
      internal_ids[tid] = extid
    }
  end
  internal_ids = idmap.invert
  allh = Hash.new {|h,k| h[k] = []}  
  IO.readlines(options[:rankfile]).each do |l|
    if /\t/ =~ l then
      l = l.split("\t")
    elsif / / =~ l then
      l = l.split(" ")
    end
    # map rank-IDs to integers
    tid = l[0].chomp
    newl=Array.new
    
    l[1..(l.size)].each do |val|
      newl << (val.chomp).to_f
    end
    allh[tid] = newl
  end

  # Choose sequences that map to rank-file, if ID is ambigous choose longest sequence
  puts "\n>> reading sequences ..." if !options[:web] 
  sequences = Hash.new
  num_reads = Hash.new
  pa_reads = false
  nrc = 0
  
  IO.readlines(options[:seqfile],">").each do |entry|
    ls = entry.split("\n").map{|x| x.chomp}
    # NEW: Remove Sequence unavailable tags in ENSEMBL sequences
    if ls[1..(ls.size)].join('') == "Sequence unavailable>" then
      next
    end
#    print "--------- sharp -------------"
    tmpl = ls[0].split("#")
    if tmpl.size > 1 then
      pa_reads = true
      ls[0] = tmpl[0]
      tmp_num_reads =  tmpl[1].to_i
      #print tmp_num_reads,"\n"
    else
      tmpl = ls[0].split("|")
      if tmpl.size > 1 then
        # 1st ID Genome
        #ls[0] = tmpl[0]
        # 2nd ID Transcript
        ls[0] = tmpl[1]
      end
    end
    # if allh.has_key?(tid)
    tid = idmap[ls[0].chomp]
    if tid then
#      print tmp_num_reads,"\n"
      if sequences.has_key?(tid) then
        newseq = ls[1..(ls.size)].join('').upcase.gsub('U','T').gsub('>','')
        if pa_reads then
#          print tmp_num_reads,"\n"
          if tmp_num_reads > num_reads[idmap[tid]].to_i then
            sequences[tid] = newseq
            nrc = nrc + 1
          elsif tmp_num_reads == num_reads[tid] && newseq.size > sequences[tid].size then
            sequences[tid] = newseq
          end
        else
        # Choose longest
          if newseq.size > sequences[tid].size then
            sequences[tid] = newseq
          end
        end
      else
        sequences[tid] = ls[1..(ls.size)].join('').upcase.gsub('U','T').gsub('>','') # last field is ">"
        num_reads[tid] = tmp_num_reads if pa_reads
        if sequences[tid].size <= 0 then
          sequences.delete(tid)
        end
      end
    end
  end
  #print "nrc ",nrc,"\n"
#  ccc = 0
#  ddd=0
#  num_reads.values.each{|v|  if v != 0 then ccc=ccc+1 else ddd=ddd+1 end}
#  print ccc," " ,ddd,"\n"
else
  # read in mirbase seed family
  word_annotation = Hash.new("") 
  if !options[:noAnn] then
    IO.readlines(annofile).each{|l| word_annotation[l.split("\t")[0]] = l.split("\t")[1]}
  end

  idmap = Hash.new
  internal_ids = Hash.new
  
  IO.readlines(tidfile).each do |l|
    tid = (l.split(" ")[0])
    l.split(" ")[1].split(",").each{|extid| 
      idmap[extid] = tid
      internal_ids[tid] = extid
    }
  end
  
  # ID to and from toFromID
  # ET - ET ; EG - EG ; EG - UC ; RS - 
  # 
  toFromID = Array.new(2)
  seqmap = Hash.new
  # read optional sequences
  idcounter = 0
  if options[:seqfile]
    puts "\n>> reading sequences ..." #if !options[:web] 
    sequences = Hash.new
  #  num_reads = Hash.new if options[:pa_seq]
    seqLengths = Hash.new
    allSequences = Hash.new
    segIDFile = ""
    firstSeq = 0
    IO.readlines(options[:seqfile],">").each do |entry|
      ls = entry.split("\n").map{|x| x.chomp}
      # hash ensures sequence ids unique
      # NEW
      tmpl = ls[0].split("|")
      if tmpl.size > 1 then
        # 1st ID G
        # ls[0] = tmpl[0]
        # 2nd ID T
          ls[0] = tmpl[1]
      end
       # Match sequence IDs
      if firstSeq == 2 then
        if ( (spepatt == "hsa" && ls[0].match(/^ENST[0-9]+/)) || (spepatt == "mmu" && ls[0].match(/^ENSTMUS[0-9]+/)) || (spepatt == "dme" && ls[0].match(/^FBtr[0-9]+/)) ) then
          # ENSG
          toFromID[0] = "ET"
          segIDFile = resdir + "genemap_" + spepatt + "_EG_ET.tsv"
          #print "Sequence ID 1 ",ls[0][0..3]," ",ls[1]
        elsif ( (spepatt == "hsa" && ls[0].match(/^ENSG[0-9]+/)) || (spepatt == "mmu" && ls[0].match(/^ENSGMUS[0-9]+/)) || (spepatt == "dme" && ls[0].match(/^FBgn[0-9]+/)) || (spepatt == "cel" && ls[0].match(/^[a-zA-Z0-9]+\.[0-9]+/)) ) then
          # ENST
          toFromID[0] = "EG"
          segIDFile = resdir + "genemap_" + spepatt + "_EG_ET.tsv"
          IO.readlines(segIDFile).each do |l|
            ll = l.split()
            tid = ll[0].chomp
            seqmap[tid] = tid
          end
        elsif spepatt == "hsa" && ls[0].match("^IPI[0-9]+") then
          # IPI
          segIDFile = resdir + "genemap_" + spepatt + "_EG_UC.tsv"
        elsif spepatt == "hsa" && ls[0].match("^uc[A-Za-z0-9]+\.[0-9]+") then
          # UCSC    
          segIDFile = resdir + "genemap_" + spepatt + "_EG_UC.tsv"
        elsif (spepatt == "hsa" && ls[0].match(/^[ANXYZ][PCGMRTWZS]_[0-9]+/)) || (spepatt == "cel" && ls[0].match(/^[ANXYZ][PCGMRTWZS]_[0-9]+/)) then
          # RefSeq
          segIDFile = resdir + "genemap_" + spepatt + "_EG_RS.tsv"
        elsif (spepatt == "hsa" || spepatt == "mmu" || spepatt == "dme") && ls[0].match(/^[A-Za-z0-9]+/) then
          # HUGO Symbol
          segIDFile = resdir + "genemap_" + spepatt + "_EG_HG.tsv"
        end
      end
      firstSeq = firstSeq + 1
      seq = ls[1..(ls.size)].join('').upcase.gsub('U','T').gsub('>','').gsub("SEQTENCE TNAVAILABLE",'') # last field is ">"
      if allSequences[ls[0].chomp].nil? then
        aSize = 0
        else
        aSize =  allSequences[ls[0].chomp].size
      end
      if seq.length > 0 && seq.length > aSize  then
        allSequences[ls[0].chomp] = seq
      end
    end
    # read in ID conversion table and choose longest sequence that maps to an ambiguos ID
    IO.readlines(segIDFile).each do |l|
      ll = l.split()
      tid = ll[0].chomp
      if toFromID[0] == "EG" then  
        extid = ll[0].chomp
      else
        extid = ll[1].chomp
      end
      if allSequences[extid].nil? then
        lSeq = 0
      else
        lSeq = allSequences[extid].size
      end
      if seqLengths[tid].nil? || seqLengths[tid] < lSeq then
        seqLengths[tid] = lSeq
        seqmap[tid] = extid 
      end
    end
        
    redun = 0
    seqmap.each do |extid,tid|
      if allSequences[extid].nil? then
        # seqmap.delete(tid)
        redun = redun + 1
      end
    end
#    print "seqmap size  ", seqmap.size," Total number of sequences left ", allSequences.size, "\n"
  end
#
  rnkmap = Hash.new
  allRanks = Hash.new
  allh = Hash.new
  filtered = 0
  internal_ids = Hash.new
  firstSeq = 0
  ii = 0
  doublette = 0
  IO.readlines(options[:rankfile]).each do |l|
    ii = ii + 1
      ls = l.split()
#    print ls[1..ls.length],"\n"
    if ls[1].nil? then
      next
    end
  #  print ls.join(" ")," ",ii,"\n"
    if allRanks[ls[0]].nil? then
      allRanks[ls[0]] = ls[1..ls.length]
    else
      #print ls[0],ls[1],"\n"
      if allRanks[ls[0]][0] < ls[1] then
        allRanks[ls[0]] = ls[1..ls.length]
        doublette = doublette + 1
      end
    end
    #print ls[1..-1], " " , ls[0..-1] , "\n"
    # Match rank IDs
    if firstSeq == 2 then
      if  ( (spepatt == "hsa" && ls[0].match(/^ENST[0-9]+/)) || (spepatt == "mmu" && ls[0].match(/^ENSTMUS[0-9]+/)) || (spepatt == "dme" && ls[0].match(/^FBtr[0-9]+/)) ) then
        # ENSG
        toFromID[1] = "ET"
        segIDFile = resdir + "genemap_" + spepatt + "_EG_ET.tsv"
        #    print "Sequence ID 1 ",ls[0][0..3]," ",ls[1]
      elsif ( (spepatt == "hsa" && ls[0].match(/^ENSG[0-9]+/)) || (spepatt == "mmu" && ls[0].match(/^ENSGMUS[0-9]+/)) || (spepatt == "dme" && ls[0].match(/^FBgn[0-9]+/)) || (spepatt == "cel" && ls[0].match(/^[a-zA-Z0-9]+\.?[0-9]+/)) ) then
        # ENST
        toFromID[1] = "EG"
        #   print "Sequence ID 2 ",ls[0][0..3]," ",ls[1]
        segIDFile = resdir + "genemap_" + spepatt + "_EG_ET.tsv"
        IO.readlines(segIDFile).each do |l|
          ll = l.split()
          tid = ll[0]
          rnkmap[tid] = tid
        end
      elsif spepatt == "hsa" && ls[0].match(/^IPI[0-9]+/) then
        # IPI
        segIDFile = resdir + "genemap_" + spepatt + "_EG_IP.tsv"
      elsif spepatt == "hsa" && ls[0].match(/^uc_[0-9.]+/) then
        # UCSC
        segIDFile = resdir + "genemap_" + spepatt + "_EG_UC.tsv"
      elsif spepatt == "hsa" && ls[0].match(/^[ANXYZ][PCGMRTWZS]_[0-9]+/) then
        # RefSeq
        segIDFile = resdir + "genemap_" + spepatt + "_EG_RS.tsv"
      elsif (spepatt == "hsa" || spepatt == "mmu" || spepatt == "dme") && ls[0].match(/^[A-Za-z0-9]+/) then
        # HUGO Symbol
        segIDFile = resdir + "genemap_" + spepatt + "_EG_HG.tsv"
      end
      if toFromID[1] != "EG" then
        IO.readlines(segIDFile).each do |lt|
          ll = lt.split()
          tid = ll[0]
          extid = ll[1]
          rnkmap[extid] = tid.to_s
        end
      end
    end
    firstSeq = firstSeq + 1
  end

#  print "rnkmap size " ,rnkmap.size," allRanks size " ,allRanks.size," Doublettes ",doublette, "\n"
  # go through sequences and map transcript ID to longest sequences
  if options[:gene_set] then
    internal_ids = Hash.new
  end
  iid = 0
  failed = 0
  unmapp = 0
  sequences = Hash.new
  allRanks.each do |extid,vals|
    # subset = Array.new
    longestSeq = ""
    longestSeqTid = ""
    currentSeq = ""
    if rnkmap[extid].nil? then
      invFile.print "No_Mapping" + "\t" + extid.to_s,"\n"
      unmapp = unmapp + 1
      currentSeq = ""
    else
      tid = rnkmap[extid]
      if allSequences[seqmap[tid]].nil? then
        invFile.print "No_Sequence" + "\t" + seqmap[tid].to_s,"\n"
        currentSeq = ""
        unmapp = unmapp + 1
      else
        currentSeq = allSequences[seqmap[tid]]
      end
    end
    if currentSeq != "" then
      if options[:gene_set] then
        internal_ids[tid] = iid
      end
      sequences[iid.to_s] = currentSeq

      allh[iid.to_s] = Array.new
      vals.each_with_index do |val,i| 
        allh[iid.to_s] << val.to_f
      end
    else
      failed = failed + 1
    end
    iid = iid + 1
    #print "int id ",iid
  end

#print "Invalid sequences ", failed-unmapp," Unmappable " ,unmapp ,"\n"
filtered = failed
end

nwords = 0
options[:wordsize].each{|ws| 
  nwords = nwords + 4**ws
}

puts "\n>> Mapping and filtering IDs ..." if !options[:web]

all = []

# filter unknown sequences
sequences.keys.each{|id| 
  if !allh.key?(id) then
    sequences.delete(id)
  end
}

allh.keys.each{|id| 
  if !sequences.key?(id) then
    allh.delete(id)
    invFile.print "No_Sequence" + "\t" + internal_ids[id].to_s,"\n"
    filtered = filtered + 1
  end
} 

# we currently mean-collapse ids, we could allow mean/min/max collapsing ...
all = allh.to_a.map{|tid,values| [tid.to_s,values.to_statarray.mean]}
allh = nil

if !options[:web] then
  puts "removed #{filtered} invalid IDs" if filtered > 0
else
  puts "results invalid #{filtered}"
end

# ID information is necesary if we want to print out target sites
if options[:gene_set] then
  idmap = internal_ids.invert
end
invFile.close

all.sort!{|a,b| a[1] <=> b[1]} 
if !options[:cooccur].nil? then
  cooccur = false
  coset = ""
  coscore = []             
end
secondPass = false
while(pass)
  
  pass = false
  if options[:cooccur] && !options[:one_plot_words].empty? then
    tmpwords = options[:one_plot_words]
    options[:one_plot_words] = []
  end
  if cooccur then
    copass = true
    options[:one_plot_words] = tmpwords if !tmpwords.nil?
  end
    ###
    ### Word enumeration (optional)
    ###
  
    wordscores = []
    
    # Freeing Memory
    JavaFreeMem.freeMem()
  
  if sequences
    puts "\n>> Enumerating words in sequences"  if !options[:web] && !secondPass
    nSeqs = sequences.size
    nGenes = nSeqs
    outFileArray = Array.new(threads)
    outFileNArray = Array.new(threads)
    #NEW check if there were IDs that mapped
    raise 'No IDs could be mapped from rank-file to sequence file. Either you\'re using an unsupported type of ID or your IDs are mapable from rank-file to sequence file, but you forgot the --custom_ids option' unless sequences.size != 0    
    (0..(threads - 1)).each{ |i|
      outFileNArray[i] = outDir + "seq_" + ranStr + i.to_s + ".tmp"
      if File::exists?(outFileNArray[i]) then
        File.delete(outFileNArray[i])
      end
      outFileArray[i] = File.open(outFileNArray[i],"w")
    }
    
    
#    print all.join(" ")

    kkk = 0
    lines = (2*nGenes)/threads
    if lines % 2 != 0 then
      lines = lines - 1
    end
    lines = lines / 2
    iit = 0
    outFileArray.each{ |outf|
      if kkk < (threads-1) then
        (0..lines).each{ |i|
          x = all[iit]
  #        print x[0].to_s," ",sequences.keys[2]," ",iit," ",kkk,"\n"
          if !sequences[x[0].to_s].nil? then
            outf.puts ">#{x[0]}"
            outf.puts sequences[x[0].to_s]
          else
            print "Break ",sequences[x[0].to_s]," ",x[0]," ", sequences[1.to_s],"\n"
            break
          end  
          iit = iit + 1
        }
        outf.close
        kkk = kkk + 1
      else
        x = all[iit]
        while !sequences[x[0].to_s].nil? do
          outf.puts ">#{x[0]}"
          outf.puts sequences[x[0].to_s]
          iit = iit + 1
          if ! all[iit].nil? then
            x = all[iit]
          else
            break
          end 
        end
        outf.close
      end
    }
    if !secondPass then
      
      semaphore4 = Mutex.new    
      options[:wordsize].each{ |ws|
        if varbg then
          bg = ws - 2
        end 
        if (ws < 10) then
          inMem = 1000
          #elsif (ws == 11)
          #  inMem = 400
        else
          inMem = 50
        end
        ha = Array(threads)
        pbar = ProgressBar.new("Data analysis",threads)  if !options[:web]
        threadA = Array.new(threads)
        bE = Array.new(threads)
          
        (0..( threads-1)).each do |ij|
          threadA[ij] = Thread.new{
            
            bE[ij] = BinomEvaluatorMult.new(ws,bg,outFileNArray[ij],ij)
            ha[ij] = bE[ij].getHashArray(options[:flatBgModel])
            ha[ij].each{ |key,value|
              semaphore4.synchronize{
                if newHash[calc_patt(key,ws)].nil? then
                  newHash[calc_patt(key,ws)] = value.to_a
                else
                  newHash[calc_patt(key,ws)] = newHash[calc_patt(key,ws)] + value.to_a
                end
              }
            } 
            pbar.inc  if !options[:web]
          }
        end
        (0..(threads-1)).each{|ij|
          threadA[ij].join()
        }
        pbar.finish  if !options[:web]
      }
    end
  end

  ###
  ### Generate list ranking
  ###
  (0..(threads - 1)).each{ |i|
    File.delete(outFileNArray[i])
  }
  analyze = []
  if options[:rank_split_median]
    # we should perhaps use an :inverse option,
    # reversing the two pos and neg lists
    med = all.map{|x| x[1]}.to_statarray.median
    pos_set = all.select{|x| x[1] > med}.sort{|a,b| b[1] <=> a[1]}
    neg_set = all.select{|x| x[1] <= med}.sort{|a,b| a[1] <=> b[1]}
    analyze = [[pos_set,'positive'],[neg_set,'negative']]
  elsif options[:rank_all] # do not split positive and negative range
    pos_set = all.sort{|a,b| b[1] <=> a[1]}
    neg_set = all.sort{|a,b| a[1] <=> b[1]}
    analyze = [[pos_set,'positive'],[neg_set,'negative']]
  elsif options[:rank_abs] # rank by absolute values
    pos_set = all.map{|x| [x[0],x[1].abs]}.sort{|a,b| b[1] <=> a[1]}
    neg_set = pos_set.reverse 
    analyze = [[pos_set,'positive'],[neg_set,'negative']]
 elsif options[:dist]
    tmparr = Array(all.length)
    all.each_with_index{|it,id| tmparr[id] = it[1]}
    sArr =  tmparr.to_statarray
    allstd =  sArr.stddev
    allMean =  sArr.mean
    pos_set = all.sort{|a,b| b[1] <=> a[1]}
    pos_set = pos_set.select{|item| (item[1]-allMean)/allstd > 0.5}
    neg_set = all.sort{|a,b| a[1] <=> b[1]}
    neg_set = neg_set.select{|item| (item[1]-allMean)/allstd < -0.5}
    analyze = [[pos_set,'positive'],[neg_set,'negative']]
  elsif options[:lead]
    pos_set = all.sort{|a,b| b[1] <=> a[1]}
    pos_set = pos_set[0..(options[:lead] - 1)]
    neg_set = all.sort{|a,b| a[1] <=> b[1]}
    neg_set = neg_set[0..(options[:lead] - 1)]
    analyze = [[pos_set,'positive'],[neg_set,'negative']]
  elsif options[:inv_lead]
    pos_set = all.sort{|a,b| b[1] <=> a[1]}
    setlength = pos_set.length
    pos_set = pos_set[(setlength-(options[:inv_lead]))..-1]
    neg_set = all.sort{|a,b| a[1] <=> b[1]}
    setlength = neg_set.length
    neg_set = neg_set[(setlength-(options[:inv_lead]))..-1]
    analyze = [[pos_set,'positive'],[neg_set,'negative']]
  elsif gene_intv
    pos_set = all[gene_intv[0]..gene_intv[1]].sort{|a,b| b[1] <=> a[1]}
    neg_set = all[gene_intv[0]..gene_intv[1]].sort{|a,b| a[1] <=> b[1]}
    analyze = [[pos_set,'positive'],[neg_set,'negative']]
  else
    pos_set = all.select{|x| x[1] > 0}.sort{|a,b| b[1] <=> a[1]}
    neg_set = all.select{|x| x[1] < 0}.sort{|a,b| a[1] <=> b[1]}
    analyze = [[pos_set,'positive'],[neg_set,'negative']]
  end

  # inverse lists
  analyze.map!{|set,nm| [set.reverse,nm+".inv"]} if options[:rank_inverse]
  
  # split sequence set when --split option is given
  if options[:split_words]
    seqs_with_words = Hash.new
    
    options[:split_words].each do |split_word|
      begin
        IO.readlines(prankdir + split_word.upcase + ".rnk").each do |x|
          l = x.split("\t")
          seqs_with_words[l[0]] = 1 if l[1].to_i > 0
        end
      rescue
        warn "could not split sequences on word #{split_word}: " + $!
      end
    end
    
    analyze_split = []
    analyze.each do |set,nm|
      analyze_split += set.partition{|x| seqs_with_words.key?(x[0])}.zip([nm+".split+"+options[:split_words].join(","),nm+".split-"+options[:split_words].join(",")])
    end
    analyze = analyze_split
  end
  
  ###
  ### Correlation analysis
  ###
  
  puts "\n>> Analyzing sequence sets: " + analyze.map{|x| x[1]}.join(", ") if !options[:web] && !secondPass
  # max number of decimal in pval regulization
  #reg = 1/(5001)
  lengths = Array(1)
  reg = 0.000001
  numberOfGenes = Array(2)
  analyze.each do |set,nm|
    plotFirst = true
    perms = []
   # options[:permutations].times{|i| perms << (0..set.size-1).to_a.shuffle}
    ngenes = set.size
    puts "\n>> Analyzing #{nm} set ...\nnumber of genes: #{ngenes}" if !options[:web] && !secondPass
    puts("dataset type #{nm}") if options[:web]  && !secondPass
    puts("dataset ngenes #{ngenes}") if options[:web]  && !secondPass
    if nm == "positive" then
      numberOfGenes[0] = ngenes
    else
      numberOfGenes[1] = ngenes
    end
    next if ngenes == 0
    report = []
    pfdrz = []
    
    franks = Hash.new # tid => index in set
    set.each_with_index{|x,i|
      franks[x[0]] = i
    }
    
    if options[:gene_set] then
      inv_franks = franks.invert
      lim = 1000
      word_seq_file = File.new(outDir + rankfilename + "_" + "ordered_" + nm + ".fa","w")
      (0..(ngenes-1)).each do |it|
        k = internal_ids[inv_franks[it]]
        word_seq_file.print ">" + k.to_s + "\n"
        str1 = sequences[idmap[k]].downcase
        str_idx = 0
        str_l = str1.length 
        while str_idx < str_l do
          plus = 49
          if (str_idx + plus) > str_l then
            plus = str_l - str_idx
          end
          word_seq_file.print str1[str_idx..(str_idx + plus)], "\n"
          str_idx = str_idx + 50
        end
      end
      word_seq_file.close
      print "\n",outDir,"\n"
    end
    
    pbar = ProgressBar.new("progress",(nwords).to_i)  if !options[:web]
    numpos=0
    widx = 0

    if options[:gene_set] then
      if (nm[/pos/]).nil? then
        gene_set_file = File.new(outDir + rankfilename + ".neg.genesets.dat","w")
      else
        gene_set_file = File.new(outDir + rankfilename + ".pos.genesets.dat","w") 
      end
    end

    if (nm[/pos/]).nil? && (!options[:plot_neg_words].empty? || !options[:plot_words].empty? || !options[:one_plot_words].empty?) then
      negplotfile = File.new(outDir + rankfilename + ".neg.plotfile.dat","w") 
      lll = (0..(ngenes + 1)).to_a.join(" ")
      negplotfile.print lll, "\n\n"
    end
    if (nm[/neg/]).nil? && (!options[:plot_pos_words].empty? || !options[:plot_words].empty? || !options[:one_plot_words].empty?) then
      posplotfile = File.new(outDir + rankfilename + ".pos.plotfile.dat","w") 
      lll = (0..(ngenes + 1)).to_a.join(" ")
      posplotfile.print lll, "\n\n"
    end
 
    semaphore = Mutex.new
    semaphore2 = Mutex.new
    semaphore3 = Mutex.new
    
    threadA = Array.new(threads)
    widBase = 0
    
    # Freeing Memory
    JavaFreeMem.freeMem()
    

    # Do word correlation analysis
    options[:wordsize].each do |ws|
      # Threaded each loop, wordIDs are distributed on all threads
      (0..( threads-1)).each do |ij|
        wordid2 = ij-threads
        
        threadA[ij] = Thread.new{
          while wordid2 < ((4**ws) - threads) do
            wordid2 = wordid2 + threads
            wid2 = widBase + wordid2
            word2 = calc_patt(wordid2,ws).downcase 
            if options[:gene_set] then 
              if !options[:gene_set].include?(word2) then
                next
              end
            end
            pbar.inc  if !options[:web]
            #only process annotated words
            next if options[:onlyanno] and not word_annotation.key?(word)
            if !options[:plot_words].empty?  || !options[:one_plot_words].empty? then
              if !options[:plot_words].include?(word2) && !options[:one_plot_words].include?(word2) then
                next
              end
            end
            
            if (nm[/pos/]).nil? && !options[:plot_neg_words].empty? then
              next if !options[:plot_neg_words].include?(word2)
            end
            if (nm[/neg/]).nil? && !options[:plot_pos_words].empty? then
              next if !options[:plot_pos_words].include?(word2)
             
            end

            # Get score distribution
            score = Array.new(ngenes,0)
            currscore = Array.new(ngenes,0)
            sdist = Array
            # synchronized access to newHash using semaphores
            semaphore.synchronize{
              sdist = newHash[word2.upcase]
              if !sdist.nil? then
                sdist.each_with_index{|it,i|
                  tmp2 = it.to_f
                  idx2 = tmp2.to_i
                  pval2 = tmp2 - idx2.to_i
                  idx2 = idx2.to_s
              
                  # Regulization of very small pvals and calculation of log-scores
                  if !franks[idx2].nil? && ngenes > franks[idx2] then
                    score[franks[idx2]] = -Math.log(pval2 + reg)
                  end
                }
              else
                next
              end
            }

            # Check if words cooccur
            if !options[:cooccur].nil? && copass then
              if !options[:cooccur].include?(word2) then
                #print coset, nm,"\n"
                if coset != nm then
                  #print "Reverse"
                  coscore.each_with_index{|x,i|
                    currscore[ngenes-i-1] = x                    
                  }
                  #currscore.reverse
                else
                  coscore.each_with_index{|x,i|
                    currscore[i] = x                    
                  }
#                  currscore = coscore
                end
                #print "currscore ",currscore,"\n"
                score.each_with_index{ |x,idx|
                  if currscore[idx] > 0 && score[idx] > 0 then
                    score[idx] = score[idx] + currscore[idx]
                  else
                    score[idx] = 0
                  end
                }
              else
                next
              end
            end
            
            # Calculate the observed running sum and descriptive statistics
            smean = score.to_statarray.mean
            maxrs = 0
            contrib_genes = 0
            leading_edge = 0
            # running sum
            rs = 0 
            score.each_with_index do |x,i|
              if x > 0 then
                contrib_genes = contrib_genes + 1
                rs += (x - smean)
                if rs.abs > maxrs.abs
                  maxrs = rs
                  leading_edge = i + 1
                end
                
              else
                rs += -smean
              end
            end
            
            
            next if maxrs <= 0 || contrib_genes <= 5
            numpos=numpos+1               
            
            # Get occurrence data of word of interest 
            if !options[:cooccur].nil? && !copass then
              if !cooccur then
                if options[:cooccur].include?(word2) then
                  #print copass, cooccur, "\n"
                  print word2,maxrs,"/n"
                  coset = nm
                  coscore = score
                  cooccur = true
                  pass = true
                else
                  next
                end
              else
                next
              end
            end
            

            if options[:gene_set] then
              if options[:gene_set].include?(word2) then
                word_seq_file = File.new(outDir + word2 + "." + nm + ".fa","w")
                gene_set_file.print("\n" + word2 + " ")
                (0..leading_edge).each{ |ii2|
                  if score[ii2] > 0 then
                    gene_set_file.print internal_ids[inv_franks[ii2]] + " "
                    word_seq_file.print ">",internal_ids[inv_franks[ii2]], "\n"
                    k = all[inv_franks[ii2].to_s]
                    word_seq_file.print sequences[k[0].to_s].downcase,"\n"
                  end
                }
               word_seq_file.close
              else
                next
              end
            end

            # we are only interested in pos. maxrs scores,
            # because we currently analyze up/down regulated seperately
       
            # Simons model 1: moment estimates of Brownian Bridge absolute max distribution (Kolmogorov Distribution)
            t = ngenes
            ps = PermutationStats.new()
            sstd = score.to_statarray.stddev
            pmean = ps.calcMean(sstd,t)
            pstd = Math.sqrt(ps.calcVar(sstd,t,pmean))    
            pval = ps.pval(maxrs,sstd, t)
            zsc = (maxrs-pmean)/pstd
            #                print "base Pmean ", pmean/sstd,"base Pstd ",pstd/sstd,"\n"
            # Simons model 2 max distribution
            #rayleighSigma = calcSig(score)
            #            pmean = (1/2.0) * sstd * Math.sqrt(t * Math::PI)
            #            pstd = Math.sqrt(sstd**2 * t - pmean**2)
            #            pval = Math.exp(-(maxrs**2/(t*sstd**2)))
            
            #            if (nm[/neg/]).nil? then
            # positive set
            #            else
            # negative set
            #              maxrs = -maxrs
            #            end
            #            zsc = (maxrs-pmean)/pstd
            #            zsc = maxrs
            #print  "Maxrs ", maxrs, " sstd " , sstd, " pmean " ,pmean, " pstd " , pstd, " pval " , pval.to_e, " zsc " ,zsc, "\n"
            
            # Plotting normalized running sum landscape plots and synchronized access to report
            semaphore2.synchronize{
              report << [wid2, zsc,pval,nil,leading_edge]
            }
            
            if (nm[/pos/]).nil? && (!options[:plot_neg_words].empty? || !options[:plot_words].empty? || !options[:one_plot_words].empty?) then
              score = score.reverse
              arCumsum = Array.new(score.size + 1)
              sum = -pmean/pstd
              arCumsum[0] = sum
              score.each_with_index{ |x,idx|
                sum += (-x+smean)              
                arCumsum[idx+1] = (sum-pmean)/pstd
              }
              semaphore3.synchronize{
                lll = word2 + " " + zsc.to_s + " " + arCumsum.map{|x| x.to_e(2)}.join(" ")
                negplotfile.print lll, "\n\n"
              }
            end
            if (nm[/neg/]).nil? && (!options[:plot_pos_words].empty? || !options[:plot_words].empty? || !options[:one_plot_words].empty?) then
              arCumsum = Array.new(score.size + 1)
              sum = -pmean/pstd
              arCumsum[0] = sum
              score.each_with_index{ |x,idx|
                sum += (x-smean)              
                arCumsum[idx+1] = (sum-pmean)/pstd
              }
              semaphore3.synchronize{
                lll = word2 + " " + zsc.to_s + " " + arCumsum.map{|x| x.to_e(2)}.join(" ")
                posplotfile.print lll, "\n\n"
              }
            end     
          end
        }
      end
      (0..(threads-1)).each{|ij|
        threadA[ij].join()
      }
      widBase = widBase + 4**ws
    end # wordsize
    pbar.finish  if !options[:web]
    
    ###
    ### FDR
    ###
    
    # Freeing Memory
    JavaFreeMem.freeMem()
    
    # Calculation of pval, FDRs and preparation of summary output
#    puts "fdr calculation ..." if !options[:web]  && !secondPass
    fdrrank=[]

    fdrrank = pfdrz.map{|x| [x,nil]} if options[:permutations] != 0 # [zscore,word_report_index] 
    report = report.sort_by{|x| x[0]}
    report.each_with_index{|x,idx| fdrrank << [x[1],idx]}
    fdrrank.sort!{|xxx,yyy| xxx[0] <=> yyy[0]}
    mmm=report.length
    report = report.sort_by{|x| x[2]}

    report.each_with_index{ |x,kkk| x[3] = (x[2]*(mmm)/(kkk+1))}
        
    cutoff_fdr = [0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.5] if !options[:web]  && !secondPass
    puts "" if !options[:web]  && !secondPass
    puts (["fdr <="] + cutoff_fdr.map{|x| x.to_s(3)} + ["total"]).join("\t") if !options[:web]  && !secondPass
    puts (["count"] + cutoff_fdr.map{|x| report.select{|y| y[3] <= x}.size} + [report.size]).join("\t") if !options[:web] && !secondPass

    ###
    ### Output summarization
    ###
    wss = options[:wordsize]
        
    report = report.sort_by{|x| x[1]}.reverse
    if options[:web]  && !secondPass then
      puts "\nTop #{output_top} words"
      puts ['rank','word','z-score','p-value','fdr','ledge','annotation'].map{|x| sprintf("%-10s",x)}.join('')
      cutoff_fdr = [0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.5]
      puts((cutoff_fdr.map{|x| x.to_s(3)} + ["total"]).join(","))
      puts((["dataset fdr_count",(cutoff_fdr.map{|x| report.select{|y| y[3] <= x}.size} + [report.size]).join(",")]).join(" "))
      
      report[0..(300-1)].each_with_index do |r,i|
        wd = calc_patt2(r[0],wss).downcase
        s = [wd,r[1].to_s(2),r[2].to_e(2),r[3].to_e(2),r[4].to_s,word_annotation[wd]]
        puts "word " + s.map{|x| sprintf("%-1s",x)}.join(';')
      end
    end

    if options[:genplot] > 0 then
      report[0..(options[:genplot]-1)].each_with_index do |r,i|
        if (nm[/pos/]).nil? then       
          negTopWords[i] = calc_patt2(r[0],wss).downcase
        else
          posTopWords[i] = calc_patt2(r[0],wss).downcase
        end
      end
    end
    if !options[:web]  && !secondPass then
      puts "\nTop #{output_top} words"
      puts ['rank','word','z-score','p-value','fdr','ledge','annotation'].map{|x| sprintf("%-15s",x)}.join(' ')
      report[0,output_top].each_with_index do |r,i|
        wd = calc_patt2(r[0],wss).downcase
        s = [i+1,wd.chomp,r[1].to_s(2),r[2].to_e(2),r[3].to_e(2),r[4].to_s,word_annotation[wd].chomp]
        puts s.map{|x| sprintf("%-15s",x)}.join(' ')
      end
    end

    if options[:report_words]
      puts "......"
      report.each_with_index do |r,i|
        wd = calc_patt2(r[0],wss).downcase
        if options[:report_words].include?(wd) # and i > output_top
          s = [i+1,wd.chomp,r[1].to_s(2),r[2].to_e(2),r[3].to_e(2),r[4].to_s,word_annotation[wd].chomp]
          puts s.map{|x| sprintf("%-15s",x)}.join(' ')
        end
      end
    end
    
    if options[:dump] >= 0 && !secondPass then
      fname = rankfilename + options[:wordsize].join("") + ".#{nm}." + options[:dump].to_s 
      of = File.new(outDir + fname,"w")
      of.puts ['rank','word','z-score','p-value','fdr','ledge','annotation'].map{|x| sprintf("%s",x)}.join(' ')
      puts "dumping top #{options[:dump]} words in file: #{fname}" #if !options[:web] 
      report[0..(options[:dump]-1)].each_with_index do |r,i|
        wd = calc_patt2(r[0],wss).downcase
        s = [i+1,wd.chomp,r[1].to_s(2),r[2].to_e(2),r[3].to_e(2),r[4].to_s,word_annotation[wd].chomp]
        of.puts s.map{|x| sprintf("%s",x)}.join(' ')
      end
      of.close
    end
    if options[:gene_set] then
      gene_set_file.close
    end
    
    if (nm[/pos/]).nil? && (!options[:plot_neg_words].empty? || !options[:plot_words].empty? || !options[:one_plot_words].empty?) then
      negplotfile.close
    end
    if (nm[/neg/]).nil? && (!options[:plot_pos_words].empty? || !options[:plot_words].empty? || !options[:one_plot_words].empty?) then
      posplotfile.close
    end
  end
  if !secondPass then
    # Make an overview plot?
    if !options[:mkplot].nil? then
      # Has a seed site for marking on the plot not been given as option then we don't want to marks any word.
      if options[:mkplot].size != 8 then
        options[:mkplot] = "llllllll"
      end
      
      isR = system("which R > /dev/null")
      isRs = system("which Rscript > /dev/null")
      if isRs then 
        # New
        cmd = "Rscript " + basedir + "lib/mkplot.R " + rankfilename + " " + options[:mkplot] + " " + options[:wordsize].join(",") + " "  + numberOfGenes.join(",") + " " + basedir + " " + outDir + " " + options[:species] + " " + ann
      else
        if isR then 
          cmd = "R --vanilla < " + basedir + "lib/mkplot.R " + rankfilename + " " + options[:mkplot] + " " + options[:wordsize].join(",") + " "  + numberOfGenes.join(",") + " " + basedir + " " + outDir + " " + options[:species] + " " + ann
        else
          print "To produce plot Rscript or R needs to be installed\n"
        end
      end
      print "Plotting Word Cluster Plot - may take some minutes","\n" if !options[:web]
      system(cmd)# + " > /dev/null")
    end
  end
    
  if options[:genplot] > 0 && options[:plot_pos_words].empty? && options[:plot_neg_words].empty? then
    print "\n Running a second pass to save data to be plotted ", "\n\n"
     if options[:cooccur].nil? then
       secondPass  = true
       options[:plot_pos_words] = posTopWords
       options[:plot_neg_words] = negTopWords
       pass = true
    else
      if copass then
        secondPass  = true
        options[:plot_pos_words] = posTopWords
        options[:plot_neg_words] = negTopWords
        pass = true
      end
    end
  end  
 
end
if options[:genplot] > 0 && !(options[:plot_pos_words].empty?) && !(options[:plot_neg_words].empty?) then
  print "Plotting graphs - may take 0.5-5 minutes","\n"
  isR = system("which R > /dev/null")
  isRs = system("which Rscript > /dev/null")
  if options[:genplot] > 15 && options[:web] then
    overview = 15
  else
    overview = options[:genplot]
  end
  if isRs then 
    cmd = "Rscript " + basedir + "lib/plot12.R " + rankfilename +  " " + outDir +  " " + overview.to_s
  else
    if isR then 
      cmd = "R --vanilla < " + basedir + "lib/plot12.R " + rankfilename + " " + outDir + " " + overview.to_s
    else
      print "To produce plot Rscript or R needs to be installed\n"
    end
  end
  system(cmd + " 1> /dev/null 2> /dev/null")
else
  if !(options[:plot_pos_words].empty?) || !(options[:plot_neg_words].empty?) || !options[:plot_words].empty? || !options[:one_plot_words].empty? then
    print "Plotting Enrichment Profile Plots...","\n"
    isR = system("which R > /dev/null")
    isRs = system("which Rscript > /dev/null")
    if options[:genplot] > 20 && options[:web] then
      overview = 20
    else
      overview = options[:genplot]
    end
    if !options[:one_plot_words].empty? then
      if isRs then 
        cmd = "Rscript " + basedir + "lib/plotOne.R " + rankfilename + " " + outDir + " " + overview.to_s
      else
        if isR then 
          cmd = "R --vanilla < " + basedir + "lib/plotOne.R " + rankfilename + " " + outDir + " " + overview.to_s
        else
          print "To produce plot Rscript or R needs to be installed\n"
        end
      end
      print(cmd)
      system(cmd + " 1> /dev/null")
    end
    if !options[:plot_words].empty? then
      if isRs then 
        cmd = "Rscript " + basedir + "lib/plot12.R " + rankfilename + " " + outDir + " " + overview.to_s
      else
        if isR then 
          cmd = "R --vanilla < " + basedir + "lib/plot12.R " + rankfilename + " " + outDir + " " + overview.to_s
        else
          print "To produce plot Rscript or R needs to be installed\n"
        end
      end
      system(cmd + " 1> /dev/null")
    end
    
  end
end

# print directory of all output
system("mv " + outDir + "*negative.jpg " + outDir + "neg.jpg 2> /dev/null") if options[:web]
system("mv " + outDir + "*positive.jpg " + outDir + "pos.jpg 2> /dev/null") if options[:web]
print "Output files in", "\n" if !options[:web]
print outDir,"\n\n" if !options[:web]

# Web server
########################################
#table_name  field_name value
