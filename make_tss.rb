STDOUT.sync = true
require './WigReader.rb'
def error_log(message)
  e = File.open(Dir.getwd + "/WigReader_error.log","a")
  e.puts "#{__FILE__} at #{Time.now} : #{message}"
  e.close
end
error_log "New analysis"

TSS_OUTSIDE = 3000    # how far does the promoter reach away from the gene body
TSS_INSIDE  = 3000    # how far inside the genebody does the promoter reach

BED_FILE = path_to_bed_file
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
def read_bed_file(bed_file)

  hash = {} # Hash for return data
  
  # read bed file containing all gene coordinates and split them into promoter, genebody and terminationsite
  File.readlines(bed_file).each do |line|
    chr,st,en,name,_,strand,_ = line.chomp.split("\t")
    st = st.to_i
    en = en.to_i

    gene_identifier = name.split("|")[1] || name.split("|")[0]  # returns MGI symbol if possible, otherwise ENSEMBL ID

    # remove all features where the genebody would be smaller than the internal regions of promoter and TTS
    if strand == "+"
      #
      # GENE:     st                                     en
      # -----------=======================================----------------  ( - 500 bp)
      # Promoter:
      #    --------=
      # Genebody:
      #             ======================================
      # TTS:
      #                                                   --
      promoter_start = st - TSS_OUTSIDE
      promoter_end   = st + TSS_INSIDE
   
    elsif strand == "-"
      #
      # GENE:     st                                     en
      # -----------=======================================----------------  ( - 500 bp)
      # Promoter:
      #                                                  =--------
      # Genebody:
      #             =====================================
      # TTS:
      #           --
      promoter_start = en - TSS_INSIDE
      promoter_end   = en + TSS_OUTSIDE
    
    else
      # skip this line if we don't know the orientation of gene
      next
    end
    hash[gene_identifier] = {
      :chr    => chr,
      :start  => promoter_start,
      :ending => promoter_end
    }
  end
  return hash
end

coordinates = read_bed_file(BED_FILE)
puts "#{__FILE__}: Finished reading .bed file, created chromosomal location index. Starting to read .wig files"
WIG_FILES.each do |file|
  label    = file[0]
  filename = file[1]
  puts "#{__FILE__}: Reading #{label}"
  wig = WigReader.new(WIG_DIR + filename)
  wig.read
  puts "#{__FILE__}: Finished indexing. Starting analysis"

  # Open a file to write results to
  File.open(label + "_tss.txt","w") do |f|
    header = %w(symbol coordinates TSS_area TSS_peak)
    f.puts header.join("\t")

    # Calculating fpkm values for TSS
    coordinates.each do |mgi_symbol,hash|
      begin
        area = wig.fpkm(hash)
        peak = wig.profile(hash).max
      rescue => e
        error_log "#{e} => #{mgi_symbol} #{hash}"
        area = 0
        peak = 0
      end
      data = [
        mgi_symbol,
        "#{hash[:chr]}:#{hash[:start]}-#{hash[:ending]}",
        area,
        peak
      ]
      f.puts data.join("\t")
    end
  end

  # free up memory space
  wig = nil
  GC.start
end