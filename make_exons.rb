STDOUT.sync = true
require './WigReader.rb'

def error_log(message)
  e = File.open(Dir.getwd + "/WigReader_error.log","a")
  e.puts "#{__FILE__} at #{Time.now} : #{message}"
  e.close
end

error_log "New analysis"


BED_FILE = path_to_bed_file
SYMBOL_LOOKUP_TABLE = path_to_symbol_lokup_table
WIG_DIR  = path_to_wig_files

WIG_FILES = [
] # [ ["Label","filename"]

WIG_FILES.each do |file|
  f = File.open(WIG_DIR + file[1])
  f.gets
  f.close
end
puts "#{__FILE__}: Able to open and read from all #{WIG_FILES.size} files."

# reads bed file and returns a Hash => {symbol => {:chr => chr1, :start => 1001, :ending => 10013}, ... }
def read_bed_file(bed_file, lookup_table)
  if !lookup_table.nil?
    # convert NM_identifier from bed file into MGI symbols
    conversions = {}
    File.readlines(lookup_table).each do |line|
      _,_,_,name,_,_,_ = line.chomp.split("\t")
      ensembl,mgi = name.split("|")
      conversions[ensembl] = mgi || ensembl || (raise RuntimeError, "No symbol or genename found!")
    end
  end
  hash = {} # Hash for return data
  
  # read bed file containing all gene coordinates and split them into promoter, genebody and terminationsite
  File.readlines(bed_file).each do |line|
    name,exon,chr,_,st,en,_ = line.chomp.split("\t")
    st = st.to_i
    en = en.to_i

    gene_identifier = lookup_table.nil? ? name : ( conversions[name] || name )

    if hash.has_key?(gene_identifier)
      # add the current exon to Array
      hash[gene_identifier] << {
        :chr    => chr,
        :start  => st,
        :ending => en }
    else
      # make new Array and add exon
      hash[gene_identifier] = [{
        :chr    => chr,
        :start  => st,
        :ending => en }]  
    end
  end
  return hash
end

coordinates = read_bed_file(BED_FILE,SYMBOL_LOOKUP_TABLE)
puts "#{__FILE__}: Finished reading .bed file, created chromosomal location index. Starting to read .wig files"
WIG_FILES.each do |file|
  label    = file[0]
  filename = file[1]
  puts "#{__FILE__}: Reading #{label}"
  wig = WigReader.new(WIG_DIR + filename)
  wig.read
  puts "#{__FILE__}: Finished indexing. Starting analysis"

  # Open a file to write results to
  File.open(label + "_exons.txt","w") do |f|
    header = %w(symbol coordinates avg_exon_peak median_exon_peak number_of_exons all_coordinates)
    f.puts header.join("\t")

    # Calculating fpkm values for TSS
    coordinates.each do |mgi_symbol,array|
      peaks = []
      coord = []
      coord_ary = []
      array.each do |hash|
        begin
          peaks << wig.profile(hash).max
          coord << hash[:start]
          coord << hash[:ending]
          coord_ary << "#{hash[:chr]}:#{hash[:start]}-#{hash[:ending]}"
        rescue
          error_log "#{e} => #{mgi_symbol} #{hash}"
        end
      end
      len = peaks.length
      if len > 0 
        starting = coord.min
        ending   = coord.max
        total = peaks.inject(:+)
        average = total.to_f / len # to_f so we don't get an integer result
        sorted = peaks.sort
        median = len % 2 == 1 ? sorted[len/2] : (sorted[len/2 - 1] + sorted[len/2]).to_f / 2

        data = [
          mgi_symbol,
          "#{array[0][:chr]}:#{starting}-#{ending}",
          average,
          median,
          len,
          coord_ary.join("|")
        ]
        f.puts data.join("\t")
      end
    end
  end
  # free up memory space
  wig = nil
  GC.start
end