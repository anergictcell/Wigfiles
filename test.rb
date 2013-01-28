require './WigReader.rb'
require "test/unit"


class TestReading < Test::Unit::TestCase
  def setup
    @a = 196866951
    @b = 196867451
    @ary = [0.01241, 0.00839, 0.00839, 0.00839, 0.00839, 0.00839, 0.01241, 0.01677, 0.01677, 0.01677, 0.01476, 0.00839, 0.00973, 0.01476, 0.02852, 0.03892, 0.04193, 0.04193, 0.03724, 0.02751]

    @wig = WigReader.new("./test.wig")
    @wig.read
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

    
  end

  def test_profiles
    assert_equal(@ary, @wig.profile(:chr => 1, start: @a, ending: @b))

    # Removing last bin from profile
    assert_equal(@ary[0..-2], @wig.profile(:chr => 1, start: @a, ending: @b-25))
  end

  def test_errors
    assert_raise(ArgumentError) {@wig.read_with_progress("string")} # Wrong Argument
    assert_raise(ArgumentError) {@wig.fpkm(chr: 28, start: @a, ending: @b)} # Wrong chromosome name
    assert_raise(ArgumentError) {@wig.fpkm(chr: 1, start: @b, ending: @a)}  # start > ending

    assert_raise(RuntimeError) {@wig.fpkm(chr: 1, start: 100, ending: @b)}  # start < recorded data
    assert_raise(RuntimeError) {@wig.fpkm(chr: 1, start: "string", ending: @b)}  # start not Integer => Will raise Error, because "string".to_i < recorded data
    assert_raise(RuntimeError) {@wig.fpkm(chr: 1, start: 1.34, ending: @a)}  # Will raise Error, because 1.34.to_i < recorded data
  end

end