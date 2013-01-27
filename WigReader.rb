=begin

written by Jonas Marcello (rumarcello@gmail.com)
can be used/modified freely for any peaceful purpose
Aknowledgment would be nice, though


DOCUMENTATION

reads .wig files and creates a mappable Hash. The data can be accessed via chromosome coordinates.
For wig file format, please go here: http://genome.ucsc.edu/goldenPath/help/wiggle.html

==> THIS SCRIPT ONLY PARSES "variableStep" WIG FILES

USAGE:

wig = Wigreader.new(path/to/wig/file)

# read all data into the Hash at once
wig.read

# read all data into the Hash at once, showing a progressmeter
wig.read_with_progress(kind)
# kind = "chr" ==> prints the chromosome name each time a new chromosome is encountered to STDOUT
# kind = "number"|"percent"|"percentage" ==> print the percentage of data parsed to STDOUT
==> Note that read_with_progress is WAY slower than read

# return fpkm value of a specified region
wig.fpkm(:chr => "chr1", :start => 10034, :end => 123534) # => 2034.245
==> fpkm takes into account if given coordinates are only part of a data window (including both given coordinates!)
  eg: WIG FILE
  variable...chr1...
  10001 100.0
  10026 50.0    => 17/25 of the value will be used only ==> 34
  10051 100.0
  ...
  123526 100.0  => 9/25 of the value will be used only ==> 32
  123551 50.0

# return profile of given coordinates
wig.profile(:chr => "chr1", :start => 10034 :end => 123534) # => [50.0, 100.0, ... , 100.0]
==> takes only whole windows, surrounding the given coordinates


=end

class WigReader
  attr_reader :data
  def initialize(filename)
    raise ArgumentError, "File #{filename} doesn't exist" unless File.exists?(filename)
    @file = File.open(filename)
    @curr_chr = ""
    @curr_bp = 1
    @data = {}
    # data = {
    #         "chr1" => {
    #                  :step => 25, :positions => {
    #                                             10001 => 12, 10026 => 28...
    #                                             }
    #                  }
    #         }
    # {"chr1"=>{:step=>25, :positions=>{10001=>16, 10026=>10, 10051=>0}}, "chr2"=>{:step=>25, :positions=>{10001=>8, 10026=>6, 10051=>19}}}
    # @data["chr1"][:positions][10001] ==> 16
    @discarded_lines = [] 
    @skip_until_header = true # start with true to make sure that there is acutally a header in the begining
    @step = 0
  end

  def read()
    while (add)
    end
  end

  def read_with_progress(kind="chr")
    case kind
    when "chr"
      while (value = add)
        if value[0,1] == "c"  # every new chromosome returns the chromosome as a string
          puts value
        end
      end
    when "number" || "percent" || "percentage"
      require '/media/bigdisk2/projects/jonas/scripts/progressmeter.rb'
      prg = ProgressMeter.new(File.path(@file))
      while (add)
        if $.%1000 == 0 then
          print prg.to_s($.)
        end
      end
    else
      raise ArgumentError, "You didn't provide a valid measure for progress!"
    end
  end


  def add()
    if !(line = @file.gets)
      return nil
    end
    line.chomp!

    if line[0,1] == "v" #e.g. variableStep chrom=chr11 span=25
      _, chr, step = line.split(" ").collect{|a| a.split("=")[1]}
      @step = step.to_i
      
      if chr == @curr_chr
        @skip_until_header = true
        # sometimes the last few bp of a chromosome are less than the current span. Thus there will be a new header defining a new span for the next entry. Since this is usually just next to the end of a chromosome, we'll discard it.
        @discarded_lines << line
        return line
      else
        @cur_chr = chr
        @data[@cur_chr] = {:step => @step, :positions => {}}
        @curr_bp = 1
        @skip_until_header = false
      return chr
      end
    end
    

    # we are on a data line, not a "new chromosome"
    # skip lines that shouldn't be recorded, but save them in an extra variable if they are needed
    if @skip_until_header
      @discarded_lines << line
      return "Discarded : #{line}"
    end

    # check if the data line has correct format
    raise ArgumentError, "Wrong data format in line #{$.} \"#{line}\"" unless line =~ /^\d+\s\d*\.?\d+$/

    # actually reading the values and adding them to cur_chr now!
    @curr_bp = @curr_bp + @step
    pos,value = line.split(" ")
    pos = pos.to_i
    while (pos > @curr_bp)
      # Wig files don't have to list 0 values, so lines might be missing. We will give fill them with 0
      @data[@cur_chr][:positions][@curr_bp] = 0.0
      @curr_bp = @curr_bp + @step
    end
    # finally adding the value to the Hash
    @data[@cur_chr][:positions][pos] = value.to_f
    return line

  end

  def fpkm(args)
    chr = chr_known(args[:chr])
    pos_start = args[:start].to_i
    pos_end = args[:end].to_i

    #
    # check for correct entry of start and end variables
    #
    unless (pos_start.is_a? Integer) && (pos_end.is_a? Integer) && (pos_end > pos_start)
      raise ArgumentError, ":start and :end have to be Integers, and :end > :start. :start: #{pos_start}, :end: #{pos_end}"
    end

    step = @data[chr][:step]    # we will be using this value a lot, so lets get it once from the Hash

    #
    # shifting the reading window upstream until we find the start coordinate that matches our wig file frame
    #
    correct_start = shift_start(chr,pos_start)
    # calculate fraction of first window to use
    first_fraction = ((correct_start + step) - pos_start) / step.to_f

    #
    # shifting the reading window upstream until we find the end coordinate that matches our wig file frame
    #
    correct_end = shift_end(chr,pos_end)
    # calculate fraction of last window to use
    last_fraction = ((pos_end+1) - correct_end) / step.to_f

    #
    # if only a fraction of one window will be used
    #     
    if correct_start == correct_end
      fraction = (pos_end - pos_start) /  step.to_f
      return @data[chr][:positions][correct_start] * fraction
    end

    #
    # summing up fpkm values
    #
    fpkm = @data[chr][:positions][correct_start] * first_fraction
    curr_pos = correct_start + step
    while curr_pos <= (correct_end - step) do
      fpkm = fpkm + @data[chr][:positions][curr_pos]
      curr_pos += step
    end
    fpkm = fpkm + (@data[chr][:positions][correct_end] * last_fraction)

    return fpkm
  end


  def profile(args)
    #
    #
    # returns an array containing the fpkm values of each windows within the specified coordingates (surrounding the coordinates)
    # wig.profile(:chr => "chr1", :start => 1, :end => 100) # => [0, 10, 3, 24]
    #
    chr = chr_known(args[:chr])
    pos_start = shift_start(chr,args[:start].to_i)  # shift the start coordinate upstream until it matches a wig file data point
    pos_end = shift_end(chr,args[:end].to_i)  # shift the end coordinate upstream until it matches a wig file data point

    # check for correct entry of start and end variables
    unless (pos_end > pos_start)
      raise ArgumentError, "The end coordinate is smaller than start coordinate. :start: #{pos_start}, :end: #{pos_end}"
    end

    # place fpkm values in array to return
    fpkm_ary = []
    curr_pos = pos_start
    while curr_pos <= (pos_end) do
      fpkm_ary << @data[chr][:positions][curr_pos]
      curr_pos += @data[chr][:step]
    end

    return fpkm_ary
  end


private
  def shift_start(chr,correct_start)
    #
    # shifting the reading window upstream until we find the start coordinate that matches our wig file frame
    #
    until @data[chr][:positions].has_key?(correct_start)
      correct_start -= 1
      raise RuntimeError, "Start coordinate is less than 1" if correct_start < 1
    end
    return correct_start
  end

  def shift_end(chr,correct_end)
    #
    # shifting the reading window upstream until we find the end coordinate that matches our wig file frame
    #
    until @data[chr][:positions].has_key?(correct_end)
      correct_end -= 1
      raise RuntimeError, "End coordinate is less than 1" if correct_end < 1
    end
    return correct_end
  end

  
  def chr_known(chr)
    #
    # check for correct entry of :chr
    #
    if !@data.has_key?(chr)
      if @data.has_key?("chr#{chr}")
        chr = "chr#{chr}"
      else
        raise ArgumentError, "Chromosome #{chr} and chr#{chr} not known"
      end
    end
    return chr
  end

end



f = "/media/bigdisk2/sequencing/wig/jonas/temp.txt"

a = WigReader.new(f)
a.read
puts "fpkm:"
puts a.fpkm(:chr => 1, :start => 5760201, :end => 5760422)
puts "profile:"
puts a.profile(:chr => 1, :start => 5760210, :end => 5760420)


=begin
t = Time.now
a = WigReader.new(f)
a.read_with_progress()
puts Time.now-t

t = Time.now
a = WigReader.new(f)
a.read()
puts Time.now-t



puts "Position 5760401"
puts a.data["chr1"][:positions][5760401]
puts "fpkm"
puts a.fpkm(:chr => 1, :start => 5760410, :end => 5760420)
=end
