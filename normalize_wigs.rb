#! /usr/bin/ruby
# assumes folder structure like
# /media/bigdisk2/sequencing/alignment/jonas
# /media/bigdisk2/sequencing/wig/jonas


require 'pathname'
def get_username(curr_dir = nil)
    return Pathname.new(curr_dir || Dir.getwd).basename
end

def get_source_files(curr_dir = nil)    # get all wig files that have not been analyzed yet
    curr_dir ||= (Dir.getwd + "/")
    input_files = []
    Dir.foreach(curr_dir) do |filename|
        next unless filename =~ /^\w.*\.wig$/    # only get .wig files
        next if filename =~ /\.fpkm\.wig/   # files have been normalized already
        if File.exists?(curr_dir + filename[0,filename.size-3] + "fpkm.wig")    # checks if this file has been analyzed already
            puts "#{filename} already analyzed"
        else
            input_files << filename     # file has not been normalized yet
        end
    end
    return input_files
end

def get_aligned_read_number(file)
    f = File.open(file)
    done = nil
    while line = f.gets do 
        if line =~ /# reads with at least one reported alignment: (\d*)\s\(.*\)/ then
            return ($1.to_i ? $1.to_i : nil)
        end
    end
    f.close
    return nil
end

def normalize(value,reads)
    f=100000.0  # for rounding
    a = value/(reads/1000000.0) # reads per million
    return ((a*f).round() / f)
end


def normalize_wig(file,reads)
    f = File.open(file)
    out_name = file[0,file.size-3] + "fpkm.wig"
    out = File.open(out_name,"w")
    while line = f.gets do
        if line =~ /variableStep/ then
            out.puts line.chomp
            puts "Reading #{line}"
            next
        end
        if line =~ /^\d/
            position,value = line.chomp.split("\t")
            norm_value = normalize(value.to_f,reads)
            out.puts "#{position}\t#{norm_value}"
        else 
            puts "There was a line in your .wig file that didn't make sense to me:\n#{line}"
            $error_log << "Invalid line in .wig file #{file}:\n#{line}"
            f.close
            out.close
            return nil
        end
    end
    f.close
    out.close
    return out_name
end

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#        PROGRAM START

curr_dir = (ARGV[0] ? ARGV[0] : nil)
$error_log = []
time = Time.new


puts "\n\n"
input_files = get_source_files(curr_dir)            # get all .wig files in the given directory
puts "Will normalize the following #{input_files.size} files:"
puts input_files
user = get_username(curr_dir)                       # get current subfolder name, hope this is the name of the user. Use it to look for the stats file in alignment folder
puts user

input_files.each do |file|
    # get stats file from user subfolder to retrieve aligned reads
    stderr_file = "/media/bigdisk2/sequencing/alignment/#{user}/#{file[0,file.size-7]}stderr.txt"        
    if File.exists?(stderr_file) then
        reads = get_aligned_read_number(stderr_file)
        if reads then 
            puts "\nWig file: #{file}"
            puts "Number of aligned reads: #{reads}"
            output = normalize_wig(file,reads)          # calculate fpkm for each value in .wig file, writes to new file, returns filename
            if output then
                puts "Saved normalized values to: #{output}"
            else
                puts "Something went wrong during the normalization of #{file}"
                $error_log << "Error during normalization of #{file}"
            end
        else
            puts "Couldn't find the number of aligned reads in the file."
            $error_log << "Unable to get number of aligned reads in #{stderr_file}"
        end
    else
        puts "No file to get number of reads for file #{file}"
        $error_log << "Missing file with number of aligned reads: #{stderr_file} for wig file: #{file}"
    end
end
puts "\nFinished!\n"

if $error_log.size > 0 then
    File.open("error_#{time.year}-#{time.month}-#{time.day}_#{time.hour}-#{time.min}.log","w") do |f|
        $error_log.each do |error|
            f.puts error
        end
    end
end
    
        
