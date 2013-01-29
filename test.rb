require './WigReader.rb'
require "test/unit"


class TestReading < Test::Unit::TestCase
  def setup
    @a = 196866951
    @b = 196867451
    @ary = [0.01241, 0.00839, 0.00839, 0.00839, 0.00839, 0.00839, 0.01241, 0.01677, 0.01677, 0.01677, 0.01476, 0.00839, 0.00973, 0.01476, 0.02852, 0.03892, 0.04193, 0.04193, 0.03724, 0.02751]

    @first_bin_chr11 = 121565176
    @last_nt_chr16 = 98319126+24
    @last_nt_chr10 = 129993226+24

    @chr19_a = 61342301
    @chr19_b = 61342401
    @last_ary_chr19 = [0.02516, 0.03086, 0.02516, 0.01409, 0.00839]

    @wig = WigReader.new("./test.wig")
    @wig.read
  end

  def test_parsing
    assert_equal(@last_nt_chr16, @wig.data["chr16"][:termination])
    assert_equal(@first_bin_chr11, @wig.data["chr11"][:begin])

    assert_equal(@last_nt_chr10, @wig.data["chr10"][:termination])
  end

  
  def test_fpkm_values
    res = @ary.inject{|x,y| x+y}
    assert_equal(res, @wig.fpkm(:chr => 1, start: @a, ending: @b))
    
    # Using different chromosome assignment, but same chromosome
    assert_equal(res, @wig.fpkm(:chr => "chr1", start: @a, ending: @b))

    # Removing one bin at the end
    res = @ary[0..-2].inject{|x,y| x+y}
    assert_equal(res, @wig.fpkm(:chr => 1, start: @a, ending: @b-25))

    # Clipping the last bin by 5
    res = ((20/25.to_f) * @ary.last) + @ary[0..-2].inject{|x,y| x+y}
    assert_equal(res, @wig.fpkm(:chr => 1, start: @a, ending: @b-5))

    # Trimming the first bin by 10
    res = ((15/25.to_f) * @ary[0]) + @ary[1..-1].inject{|x,y| x+y}
    assert_equal(res, @wig.fpkm(:chr => 1, start: @a+10, ending: @b))

    # Pulling only one bin
    assert_equal(@ary[0], @wig.fpkm(:chr => 1, start: @a, ending: @a+25))

    # Pulling only part of one bin
    assert_equal(@ary[0]*(10/25.to_f), @wig.fpkm(:chr => 1, start: @a+5, ending: @a+15))

    # pull data downstream of recorded data
    res = @last_ary_chr19.inject{|x,y| x+y}
    assert_equal(res, @wig.fpkm(chr: 19, start: @chr19_a, ending: @chr19_b+50))
    # remove last bin of a chromosome
    res = @last_ary_chr19[0..-2].inject{|x,y| x+y}
    assert_equal(res, @wig.fpkm(chr: 19, start: @chr19_a, ending: @chr19_b))
    

  end

  def test_profiles
    assert_equal(@ary, @wig.profile(:chr => 1, start: @a, ending: @b))

    # require data downstream of recorded data
    assert_equal(@last_ary_chr19, @wig.profile(:chr => 19, start: @chr19_a, ending: @chr19_b+50))
    # same data but only to beginning of last bin
    assert_equal(@last_ary_chr19, @wig.profile(:chr => 19, start: @chr19_a, ending: @chr19_b+1))

    # Removing last bin from profile
    assert_equal(@ary[0..-2], @wig.profile(:chr => 1, start: @a, ending: @b-25))
  end

  def test_errors
    assert_raise(ArgumentError) {@wig.read_with_progress("string")} # Wrong Argument
    assert_raise(ArgumentError) {@wig.fpkm(chr: 28, start: @a, ending: @b)} # Wrong chromosome name
    assert_raise(ArgumentError) {@wig.fpkm(chr: 1, start: @b, ending: @a)}  # start > ending
    assert_raise(ArgumentError) {@wig.fpkm(chr: 19, start: @last_nt_chr10+1, ending: @last_nt_chr10+50)}

    assert_raise(RuntimeError) {@wig.fpkm(chr: 1, start: 100, ending: @b)}  # start < recorded data
    assert_raise(RuntimeError) {@wig.fpkm(chr: 1, start: "string", ending: @b)}  # start not Integer => Will raise Error, because "string".to_i < recorded data
    assert_raise(RuntimeError) {@wig.fpkm(chr: 1, start: 1.34, ending: @a)}  # Will raise Error, because 1.34.to_i < recorded data

  end

end