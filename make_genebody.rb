STDOUT.sync = true
require './WigReader.rb'
def error_log(message)
  e = File.open(Dir.getwd + "/WigReader_error.log","a")
  e.puts "#{__FILE__} at #{Time.now} : #{message}"
  e.close
end

error_log "New analysis"

# Read in bed file and split the file into TSS, genebody and TTS
TSS_OUTSIDE = 2000 # how far does the promoter reach away from the gene body
TSS_INSIDE  = 500  # how far inside the genebody does the promoter reach
TTS_OUTSIDE = 1000 # how far does the termination site extend behind the genebody
TTS_INSIDE  = 0 

BED_FILE = path_to_bed_file
WIG_DIR  = path_to_wig_folder

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
    next if st + TSS_INSIDE + TTS_INSIDE > en

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
      body_start     = promoter_end
      body_end       = en - TTS_INSIDE
      tts_start      = body_end
      tts_end        = en + TTS_OUTSIDE
   
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
      body_start     = st + TTS_INSIDE
      body_end       = promoter_start
      tts_start      = st - TTS_OUTSIDE
      tts_end        = body_start
    
    else
      # skip this line if we don't know the orientation of gene
      next
    end
    hash[gene_identifier] = {
      :tss => {
        :chr    => chr,
        :start  => promoter_start,
        :ending => promoter_end 
      },
      :genebody => {
        :chr    => chr,
        :start  => body_start,
        :ending => body_end 
      },
      :tts      => {
        :chr    => chr,
        :start  => tts_start,
        :ending => tts_end 
      },
      :original => {
        :chr    => chr,
        :start  => st,
        :ending => en
      }
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
  File.open(label + "_genebody.txt","w") do |f|
    header = %w(symbol coordinates TSS_area TSS_peak Genebody_average TTS_area TTS_peak)
    f.puts header.join("\t")

    # Calculating fpkm values for TSS
    coordinates.each do |mgi_symbol,hash|
      begin
        # Calculate peak height and peak area for TSS and TTS (sizes are constant across genes so we can use raw values)
        tss_area = wig.fpkm(hash[:tss])
        tss_peak = wig.profile(hash[:tss]).max
        tts_area = wig.fpkm(hash[:tts])
        tts_peak = wig.profile(hash[:tts]).max

        # For the genebody, calculate the peak area and divide by length (in theoretical bins) to get average signal / bin
        gb_area       = wig.fpkm(hash[:genebody])
        gb_width      = hash[:genebody][:ending] - hash[:genebody][:start]
        gb_bin_number = gb_width.to_f / 25
        gb_avg        = gb_area / gb_bin_number
      rescue => e
        error_log "#{e} => #{mgi_symbol} #{hash}"
        tss_area = 0
        tss_peak = 0
        tts_area = 0
        tts_peak = 0
        gb_avg   = 0
      end
      data = [
        mgi_symbol,
        "#{hash[:original][:chr]}:#{hash[:original][:start]}-#{hash[:original][:ending]}",
        tss_area,
        tss_peak,
        gb_avg,
        tts_area,
        tts_peak
      ]
      f.puts data.join("\t")
    end
  end

  # free up memory space
  wig = nil
  GC.start
end
