=begin

written by Jonas Marcello (rumarcello@gmail.com)
can be used/modified freely for any peaceful purpose
Aknowledgment would be nice, though


DOCUMENTATION

reads .wig files and creates a mappable Hash. The data can be accessed via chromosome coordinates.
For wig file format, please go here: http://genome.ucsc.edu/goldenPath/help/wiggle.html

==> THIS SCRIPT ONLY PARSES "variableStep" WIG FILES

USAGE:

wig = WigReader.new(path/to/wig/file)

# read all data into the Hash at once
wig.read

# read all data into the Hash at once, showing a progressmeter
wig.read_with_progress(kind)
# kind = "chr" ==> prints the chromosome name each time a new chromosome is encountered to STDOUT
# kind = "number"|"percent"|"percentage" ==> print the percentage of data parsed to STDOUT
==> Note that read_with_progress is WAY slower than read


# to retrieve bins etc, please note:
# The start coordinate will be INcluded, 
# the end coordinate will be EXcluded
# 1-10 means [1,2,3,4,5,6,7,8,9]


# return fpkm value of a specified region
wig.fpkm(:chr => "chr1", :start => 10034, :ending => 123534) # => 2034.245
==> fpkm takes into account if given coordinates are only part of a data window (excluding the end coordinate!)
  eg: WIG FILE
  variable...chr1...
  10001 100.0
  10026 50.0    => 17/25 of the value will be used only ==> 34
  10051 100.0
  ...
  123526 100.0  => 9/25 of the value will be used only ==> 32
  123551 50.0

# return profile of given coordinates
wig.profile(:chr => "chr1", :start => 10034 :ending => 123534) # => [50.0, 100.0, ... , 100.0]
==> takes only whole windows, excluding the end coordinate.
wig.profile(:chr => "chr1", :start => 10034 :ending => 123551) # => [50.0, 100.0, ... , 100.0]



=end

class WigReader
  attr_reader :data
  def initialize(filename)
    raise ArgumentError, "File #{filename} doesn't exist" unless File.exists?(filename)
    @file = File.open(filename)
    @curr_chr = ""
    @curr_bp = nil
    @data = {}
    # data = {
    #         "chr1" => {
    #                  :step => 25, 
    #                  :begin => 10001,
    #                  :termination => 10203025, 
    #                  :positions => {
    #                                10001 => 12, 10026 => 28...
    #                                 }
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
    return self
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

      if step.nil?  # TODO: Enable Wigreader to determine the span by itself by reading the next two lines
        raise ArgumentError, "No span value is given for chromosome #{chr} in wig file"
      end
      @step = step.to_i
  
      #
      # sometimes the last few bp of a chromosome are less than the current span.
      # Thus there will be a new header defining a new span for the next entry.
      # Since this is usually just next to the end of a chromosome, we'll discard it. 
      #      
      if chr == @curr_chr
        @skip_until_header = true
        @discarded_lines << line
        return line
      else

        # New chromosome
        unless @curr_bp.nil?  # In which case we are at the start of the wig file
          # Save the last position to the chr Hash as :termination
          @data[@curr_chr][:termination] = @curr_bp+24
        end
        # set instance variables for the new chromosome Hash
        @curr_chr = chr
        @data[@curr_chr] = {:step => @step, :positions => {}}
        @curr_bp = nil
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

    #
    # actually reading the values and adding them to curr_chr now!
    #
    pos,value = line.split(" ")
    pos = pos.to_i
    # values for chromosomes don't start at 1, but WAY later. Don't fill the array with 0s
    if @curr_bp.nil? 
      @curr_bp = pos
      @data[@curr_chr][:begin] = pos
    else
      @curr_bp += @step
    end

    # Wig files don't have to list 0 values, so lines might be missing. We will fill them with 0  
    while (pos > @curr_bp)
      @data[@curr_chr][:positions][@curr_bp] = 0.0
      @curr_bp = @curr_bp + @step
    end


    #
    # finally adding the value to the Hash
    #
    @data[@curr_chr][:positions][pos] = value.to_f
    return line

  end

  def fpkm(args)
    chr, pos_start, pos_end = check_coordinates(args)

    step = @data[chr][:step]    # we will be using this value a lot, so lets get it once from the Hash

    #
    # shifting the reading window upstream until we find the start coordinate that matches our wig file frame
    #
    correct_start = align_coordinates_with_bins(chr, pos_start)
    # calculate fraction of first window to use
    first_fraction = ((correct_start + step) - pos_start) / step.to_f



    #
    # shifting the reading window upstream until we find the end coordinate that matches our wig file frame
    #
    correct_end = align_coordinates_with_bins(chr, pos_end)
    # calculate fraction of last window to use
    last_fraction = (pos_end - correct_end) / step.to_f

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
    # returns an array containing the fpkm values of each windows within the specified coordinates
    # wig.profile(:chr => "chr1", :start => 1, :ending => 100) # => [0, 10, 3]
    #
    chr, pos_start, pos_end = check_coordinates(args)

    pos_start = align_coordinates_with_bins(chr, pos_start)  # shift the start coordinate upstream until it matches a wig file data point
    pos_end = align_coordinates_with_bins(chr, pos_end-1)  # shift the end coordinate upstream until it matches a wig file data point
    # Subtracting 1 from end coordinate to exclude the last given nucleotide

    # check for correct entry of start and end variables
    # TODO:
    # Currently small profiles (size < step) raise Error. Fix that 
    # FIXED, BUT HAS TO BE TESTED
    if (pos_end < pos_start)
      raise ArgumentError, "The end coordinate is smaller than start coordinate. :start: #{pos_start}, :ending: #{pos_end}"
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
  def check_coordinates(args)
    chr = chr_known(args[:chr])
    pos_start = args[:start].to_i
    pos_end = args[:ending].to_i

    unless (pos_start.is_a? Integer) && (pos_end.is_a? Integer)
      raise ArgumentError, ":start and :ending have to be Integers! :start: #{pos_start}, :ending: #{pos_end}"
    end

    unless (pos_end > pos_start)
      raise ArgumentError, ":ending > :start!  :start: #{pos_start}, :ending: #{pos_end}"
    end

    if pos_start > ( @data[chr][:termination] - 24 )
      raise ArgumentError, "The given :start coordinate is downstream of the recorded chromosome data #{pos_start}"
    end

    # given position is downstream of recorded data, so we just shift it to the last nt
    if (pos_end > @data[chr][:termination])
      pos_end = ( @data[chr][:termination] + 1)
    end

    return [ chr, pos_start, pos_end ]
  end

  def align_coordinates_with_bins(chr,pos)
    #
    # shifting the reading window upstream until we find the start coordinate that matches our wig file frame
    #
    until @data[chr][:positions].has_key?(pos)
      pos -= 1
      raise RuntimeError, "Coordinate is too small to align it with bins." if (pos < @data[chr][:begin])
    end
    return pos
  end
  
  def chr_known(chr)
    #
    # check for correct entry of chr
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