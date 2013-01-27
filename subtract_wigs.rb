#! /usr/bin/ruby
# version 1.0 
# 4/24/12
# Jonas Marcello @ Tarakhovsky lab, Rockefeller University
require '/media/bigdisk2/projects/jonas/scripts/progressmeter.rb'

STDOUT.sync = true
$error_log = []
$headers = /^variableStep chrom=(.*)\sspan=.*/      # keep this as a regex
$var_span = /span=25/                               # keep this as a regex
$print_negatives = true

time = Time.new

def usage()
    puts "\nDoes exactly what the name says: It takes 2 wig files, subtracts the values from each other and generates a new wig file"
    puts "Current settings:"
    puts "Header line: #{$headers}"
    puts "span: #{$var_span}"
    puts "Printing negative values: #{$print_negatives}\n\n"
    puts "Usage\t subtract_wigs.rb [minuend_file] [subtrahend_file] [output_file]"
    puts "Eg\t subtract_wigs.rb H3K4me3.wig input.wig H3K4me3-input.wig\n\n"
    exit
end

def split_path_file(filename)
    if filename =~ /\// then    # filename contains folders
        path = File.dirname(filename) + "/"
        if !(path =~ /^\//) then       # relative path given, make absolute
            path = Dir.getwd + "/" + path
        end
        file = File.basename(filename)
    else
        path = Dir.getwd + "/"
        file = filename
    end
    return [path,file]
end
        
def files_match(a,b)
    if a == b then
        return a[0].chomp
    else
        return nil
    end
    
end
        
def mismatch(out)
    puts "Hey, come on. You're comparing apples and oranges here. That doesn't work too well with me. Please use only wig files that have the same variableStep. And the same organism and reference genome. I know, it's cool to comapare zebrafish with humans, but I won't do that for you. Either you fix it for yourself or you're out of luck. Have a good day!"
    out.close
    exit
end

def work_on_lines(line,f)
    if line then
        if line =~ $headers then        # header line
            chr = $1
            if line =~ $var_span then   # header line for chromosome with the right bp span
                puts chr
                return [line.chomp, ""]
            else                        # all values that don't have the bp span specified will be discarded
                print "Skipping file #{File.basename(f.path)}\non line: #{line}"
                #f.gets
                line = f.gets
                while !((line =~ $headers) && (line =~ $var_span)) do        # pull new lines in the file until we reach a new chromosome header ==> if the new chromosome has an arbitrary span, we're fucked!
                    line = f.gets
                    print "."
                end
                puts "\nuntil: #{line}"
                return [line.chomp, ""]             # reached a new chromosome, back to buisness
            end
        else
            values = line.chomp.split("\s")
            return [values[0].to_i, values[1].to_f]
        end
    else
        return ["End of file", "End of file"]
    end
end

def print_row(lines,out)
    if (lines[0][1].is_a? Float) && (lines[1][1].is_a? Float) then
        f=100000.0  # for rounding
        subtraction = (((lines[0][1] - lines[1][1])*f).round() / f)
        if ($print_negatives || subtraction > 0) then 
            out.print lines[0][0]
            out.puts "\t#{subtraction}"
        end
    else
        out.puts lines[0][0] unless (lines[0][0] == "End of file")
    end
end

def align_coordinates(lines,f1,f2,out,positions)        # pulls new line from the file that lags behind 
    minimum = positions.map{|a| a.is_a?(String) ? 3e20 : a}.min
    if positions[0] == minimum then
        # minuend is smaller, print minuend value
        out.puts "#{lines[0][0]}\t#{lines[0][1]}"
        lines[0] = work_on_lines(f1.gets,f1)
    elsif positions[1] == minimum then
        # subtrahend is smaller, print 0 - subtrahend (or nothing, dependent on settings)
        if $print_negatives then
            out.puts "#{lines[1][0]}\t-#{lines[1][1]}"
        end
        lines[1] = work_on_lines(f2.gets,f2)
    else
        $error_log << "Failure to align coordinates.\nCurrent positions are #{lines[0][0]} & #{lines[1][0]}"
        out.close
        exit
    end
    return lines
end

if ARGV.size < 3 then
    puts "You have to specify the files you want to subtract from each other\n"
    usage()
end

minuend = split_path_file(ARGV[0])
subtrahend = split_path_file(ARGV[1])
dst_file = split_path_file(ARGV[2])


puts "Subtract:\n#{minuend}\n-\n#{subtrahend}\nSave to:\n#{dst_file}"
puts "Writing negative values: #{$print_negatives}"

prg = ProgressMeter.new(minuend[0] + minuend[1])
f1 = File.open(minuend[0] + minuend[1])
f2 = File.open(subtrahend[0] + subtrahend[1])

File.open(dst_file[0] + dst_file[1],"w") do |out|
    lines = []
    lines[0] = work_on_lines(f1.gets,f1)
    lines[1] = work_on_lines(f2.gets,f2)
    (files_match(lines[0],lines[1]) || mismatch(out))
    until lines.uniq == [["End of file", "End of file"]] do
        until (lines[0][0] == lines[1][0]) do                       # does this until both files are synced
            positions = lines.map{|line| line[0]}
            # debugging: puts "Current position call #{positions}"
            lines = align_coordinates(lines,f1,f2,out,positions)     # move along in one file, keep the other waiting to align
        end                                                         # gets here only if both files are synced
        print_row(lines,out)
        lines[0] = work_on_lines(f1.gets,f1)
        lines[1] = work_on_lines(f2.gets,f2)
        if $.%10000 == 0 then
          print prg.to_s($.)
        end
    end
end

if $error_log.size > 0 then
    File.new("error_#{time.year}-#{time.month}-#{time.day}_#{time.hour}-#{time.minute}.log","w") do |f|
        $error_log.each do |error|
            f.puts error
        end
    end
end

puts "Finished"
