STDOUT.sync = true
require_relative './WigReader.rb'
def error_log(message)
  e = File.open(Dir.getwd + "/WigReader_error.log","a")
  e.puts "#{__FILE__} at #{Time.now} : #{message}"
  e.close
end

error_log "New analysis"

if ARGV.size != 3 then
  puts "Usage:"
  puts "\t#{__FILE__} src.wig genes.bed result.profile"
  exit
end

# Read in bed file and split the file into TSS, genebody and TTS
TSS_OUTSIDE = 1000 # how far does the promoter reach away from the gene body
TSS_INSIDE  = 1000  # how far inside the genebody does the promoter reach
TTS_OUTSIDE = 1000 # how far does the termination site extend behind the genebody
TTS_INSIDE  = 1000
TSS_BINS = 20
TTS_BINS = 20
GENEBODY_BINS = 40

WIG_FILE = ARGV.shift
BED_FILE = ARGV.shift
OUT_FILE = ARGV.shift


def check_file(filename,create)
  if create then
    f = File.open(filename,"w")
    f.puts "Testing to see if we can write to the file"
    f.close
  end

  f = File.open(filename)
  if f.gets then
    f.close
    if (create) then
      # Deleting the contents of the newly created file to have an empty file
      f = File.open(filename,"w")
      f.close
    end

    return true
  else
    return false
  end
end
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
        :ending => en,
        :strand => strand
      }
    }    
  end
  return hash
end

unless check_file(WIG_FILE,false)
  puts "Unable to read from #{WIG_FILE}"
  exit
end

unless check_file(BED_FILE,false)
  puts "Unable to read from #{BED_FILE}"
  exit
end

unless check_file(OUT_FILE,true)
  throw "Unable to read/write from/to #{OUT_FILE}"
end

coordinates = read_bed_file(BED_FILE)
puts "#{__FILE__}: Finished reading .bed file, created chromosomal location index. Starting to read .wig files"

puts "#{__FILE__}: Reading #{WIG_FILE}"
wig = WigReader.new(WIG_FILE)
wig.read
puts "#{__FILE__}: Finished indexing #{WIG_FILE}. Starting analysis"

File.open(OUT_FILE,"w") do |f|
  header = %w(symbol coordinates length bins binlength)
  1.upto(TSS_BINS+GENEBODY_BINS+TTS_BINS) do |i|
    header << "BIN_#{i}"
  end
  f.puts header.join("\t")
  coordinates.each do |mgi_symbol,hash|
    tss_length = hash[:tss][:ending]-hash[:tss][:start]
    genebody_length = hash[:genebody][:ending]-hash[:genebody][:start]
    tts_length = hash[:tts][:ending]-hash[:tts][:start]

    # By rounding the number of bp in each bin we
    # accept that some genebodies will be cut up to 25 bp short
    # while other genebodies will reach up to 25 bp into the TTS
    # The shift of +/- 25 bp is acceptable, as we're talking about average
    # genebody length of 52,000 bp and average genebody bin size of > 1000 bp
    tss_bin_length = (tss_length.to_f / TSS_BINS).round
    genebody_bin_length = (genebody_length.to_f / GENEBODY_BINS).round
    tts_bin_length = (tts_length.to_f / TTS_BINS).round
    values = []
    begin
      if hash[:original][:strand] == "+" then
        (0...TSS_BINS).each do |i|
          values << wig.fpkm({
            :chr => hash[:tss][:chr],
            :start => hash[:tss][:start]+ i*tss_bin_length,
            :ending => hash[:tss][:start]+ i*tss_bin_length + tss_bin_length
            })/tss_bin_length
        end
        (0...GENEBODY_BINS).each do |i|
          values << wig.fpkm({
            :chr => hash[:genebody][:chr],
            :start => hash[:genebody][:start]+ i*genebody_bin_length,
            :ending => hash[:genebody][:start]+ i*genebody_bin_length + genebody_bin_length
            })/genebody_bin_length
        end
        (0...TTS_BINS).each do |i|
          values << wig.fpkm({
            :chr => hash[:tts][:chr],
            :start => hash[:tts][:start]+ i*tts_bin_length,
            :ending => hash[:tts][:start]+ i*tts_bin_length + tts_bin_length
            })/tts_bin_length
        end
      else 
        # reverse each element
        (0...TSS_BINS).each do |i|
          rev = TSS_BINS
          values << wig.fpkm({
            :chr => hash[:tss][:chr],
            :start => hash[:tss][:ending] - i*tss_bin_length - tss_bin_length,
            :ending => hash[:tss][:ending] - i*tss_bin_length
            })/tss_bin_length
        end
        (0...GENEBODY_BINS).each do |i|
          values << wig.fpkm({
            :chr => hash[:genebody][:chr],
            :start => hash[:genebody][:ending] - i*genebody_bin_length - genebody_bin_length,
            :ending => hash[:genebody][:ending] - i*genebody_bin_length
            })/genebody_bin_length
        end
        (0...TTS_BINS).each do |i|
          values << wig.fpkm({
            :chr => hash[:tts][:chr],
            :start => hash[:tts][:ending] - i*tts_bin_length - tss_bin_length,
            :ending => hash[:tts][:ending] - i*tts_bin_length
            })/tts_bin_length
        end
      end
    rescue => e
      error_log "#{e} => #{mgi_symbol} #{hash}"
      values = []
      1.upto(TSS_BINS+GENEBODY_BINS+TTS_BINS) do |i|
        values << 0
      end
    end
    data = [
      mgi_symbol,
      "#{hash[:original][:chr]}:#{hash[:original][:start]}-#{hash[:original][:ending]}",
      "#{tss_length}|#{genebody_length}|#{tss_length}",
      "#{TSS_BINS}|#{GENEBODY_BINS}|#{TTS_BINS}",
      "#{tss_bin_length}|#{genebody_bin_length}|#{tss_bin_length}",
      values.join("\t")
    ]
    f.puts data.join("\t")
  end
end






