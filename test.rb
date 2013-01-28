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
    assert_equal(0.38077, @wig.fpkm(:chr => 1, :start => @a, :end => @b))
    assert_equal(0.38077, @wig.fpkm(:chr => 1, :start => @a, :end => @b))
  end

  def test_profiles
    assert_equal(@ary, @wig.profile(:chr => 1, :start => @a, :end => @b))
    assert_equal(@ary[0..-2], @wig.profile(:chr => 1, :start => @a, :end => @b-25))
  end
end