Here you will find some Ruby scripts for dealing with .wig files  
This is helpful for ChIP-seq analysis  
Don't laugh, the code is probably horrible, I'm learning by doing
Some of the scripts are older, but I'm uploading them now without checking if they all work.  
Use at your own risk

- WigReader.rb:   
class to read .wig files into memory. Data is stored in "coordinates"=>"value" Hashes for quick access.  
Eventually the other scripts will be merged in here as well to have one class for everything.

- normalize_wigs.rb:  
script to convert raw values in .wig files to fpkm values. This script only works with a special folder strucure that is explained in the script.

- subtract_wigs.rb:  
script to subtract values from one .wig file from values of another. Good to remove background noise

- test.rb:  
Unit testing for WigReader.rb  

- test.wig:  
Test wig file  



