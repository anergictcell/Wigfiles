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
Documentation TBA

### make_exons.rb
Documentation TBA

### make_genebody.rb
Documentation TBA

### make_tss.rb
Documentation TBA

### tssProfile.rb
Documentation TBA

### normalize_wigs.rb
script to convert raw values in .wig files to fpkm values. This script only works with a special folder strucure that is explained in the script.
The calculations are not floating-number-safe and should NOT be used for publication.


### subtract_wigs.rb
script to subtract values from one .wig file from values of another. Good to remove background noise, like IgG. The data is subtracted, not divided! Output of negative values can be suppressed by setting `$print_negatives` to `false`.

### test.wig:
.wig file to be used for testing.



