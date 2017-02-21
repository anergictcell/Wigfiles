# WigReader

## A library for working with ChIP-seq or RNA-seq data in .wig file format

### Use at your own risk!
The library was created many years ago and was actually my first Ruby project. Thus it's not written too well and its definitely not the best tool for the job nowadays. However, it does work if you have the correct files. It doesn't include much error-checking, so while it checks the most common errors, there is a small change you might end up with the program crashing after a 4hour analysis..... 


### WigReader
A small class to parse and index .wig files. Can be used to retrieve custom signal profiles. 

In order to run efficiently, you need to have enough RAM. It reads the whole .wig file into memory, so you should have more RAM than the size of your .wig file.

Since it's reading the whole file into memory, it only makes sense to use if you need it to return FPKMs of many areas at once.

General usage:
```Ruby
wig = WigReader.new(path/to/wig/file)

# read all data into the Hash at once
wig.read

# return fpkm value of a specified region
wig.fpkm(:chr => "chr1", :start => 10034, :ending => 123534) # => 2034.245

# return profile of given coordinates
wig.profile(:chr => "chr1", :start => 10034 :ending => 123534) # => [50.0, 100.0, ... , 100.0]
```

If you are worried about the script crashing, you can read in the wig file with a progress bar. It will noticeable slow down the parsing though.
```Ruby
wig = WigReader.new(path/to/wig/file)
wig.read_with_progress(<type>)

# Output each new chromosome during parsing (default)
wig.read_with_progress("chr")

# Output the percentage of parsed data
wig.read_with_progress("percent")
```

Testing is done via `test.rb`.

### geneProfile.rb
Calculates the FPKM profiles of genomic areas. 
Takes each given area, divides it into equal bins, then calculates the fpkm value for each bin.

General usage:
```Bash
geneProfile.rb [source-wig-file] [areas-bed-file] [output-file]

geneProfile.rb ./H3K4me3.wig ./Protein.coding.genes.bed ./Protein.coding.gened.profiles.tsv
```
Some settings can be modified inside of the script:
`TSS_OUTSIDE` : Number of bp upstream of TSS that will be included for the promoter (Default: 1000)
`TSS_INSIDE` : Number of bp downstream of TSS (inside genebody) that will be included for the promoter (Default: 1000)
Default Promoter: TSS +/- 1 kb
`TTS_OUTSIDE` : Number of bp downstream of TTS that will be included for transcriptional termination site (Default: 1000)
`TTS_INSIDE` : Number of bp upstream of TTS (inside genebody) that will be included for transcriptional termination site (Default: 1000)
Default TTS: TTS +/ 1 kb

Based on the standard settings, the first and last 1 kb of the genebody will NOT be included in the genebody, but only in the TSS and TTS.

`TSS_BINS`, `TTS_BINS` : Defines in how many bins the promoter and TTS areas will be split for profile calculation (Default: 20)
`GENEBODY_BINS` : Defines in how many bins the genebody will be split (Default: 40)

### make_exons.rb
Documentation TBA

### make_genebody.rb
Documentation TBA

### make_tss.rb
Documentation TBA

### tssProfile.rb
Calculates the FPKM profiles of promoters (TSS)
Takes each given gene, divides the promoter area into equal bins, then calculates the fpkm value for each bin.

General usage:
```Bash
tssProfile.rb [source-wig-file] [areas-bed-file] [output-file]

tssProfile.rb ./H3K4me3.wig ./Protein.coding.genes.bed ./Protein.coding.gened.profiles.tsv
```
Some settings can be modified inside of the script:
`TSS_OUTSIDE` : Number of bp upstream of TSS that will be included for the promoter (Default: 3000)
`TSS_INSIDE` : Number of bp downstream of TSS (inside genebody) that will be included for the promoter (Default: 3000)
Default Promoter: TSS +/- 3 kb
`TSS_BINS` : Defines in how many bins the promoter and TTS areas will be split for profile calculation (Default: 80)


### normalize_wigs.rb
script to convert raw values in .wig files to fpkm values. This script only works with a special folder strucure that is explained in the script.
The calculations are not floating-number-safe. If your data is sensitive to those errors, do NOT use this tool.


### subtract_wigs.rb
script to subtract values from one .wig file from values of another. Good to remove background noise, like IgG. The data is subtracted, not divided! Output of negative values can be suppressed by setting `$print_negatives` to `false`.

### test.wig:
.wig file to be used for testing.



