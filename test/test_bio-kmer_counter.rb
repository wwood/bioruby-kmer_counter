require 'helper'
require 'tempfile'

class TestBioKmerCounter < Test::Unit::TestCase
  should 'test_lowest_lexigraphical_form' do
    assert_equal Bio::Sequence::NA.new('AA'), Bio::Sequence::NA.new('AA').lowest_lexigraphical_form
    assert_equal Bio::Sequence::NA.new('AA'), Bio::Sequence::NA.new('TT').lowest_lexigraphical_form
    assert_equal Bio::Sequence::NA.new('AG'), Bio::Sequence::NA.new('CT').lowest_lexigraphical_form
  end
  
  should 'test_empty_full_kmer_hash' do
    answer = {}; %w(A C G T).each{|k| answer[k] = 0}
    assert_equal answer, Bio::Sequence::Kmer.empty_full_kmer_hash(1)
  end
  
  should 'test merge down' do
    answer = {}; %w(A C).each{|k| answer[k] = 0}
    full = Bio::Sequence::Kmer.empty_full_kmer_hash(1)
    assert_equal answer, Bio::Sequence::Kmer.merge_down_to_lowest_lexigraphical_form(full)
    full = Bio::Sequence::Kmer.empty_full_kmer_hash #defaults to kmer hash length 4
    assert_equal 136, Bio::Sequence::Kmer.merge_down_to_lowest_lexigraphical_form(full).length
  end
  
  def script_path
    File.join(File.dirname(__FILE__),'..','bin','kmer_counter.rb')
  end
  
  should 'test_running1' do
    Tempfile.open('one') do |tempfile|
      tempfile.puts '>one'
      tempfile.puts 'ACAGT'
      tempfile.close

      assert_equal "ID\tA\tC\none_0\t0.6\t0.4\n", `#{script_path} -w 5 -k 1 #{tempfile.path}`
    end
  end
  
  should 'not whack out when there isnt any sequence to count' do
    Tempfile.open('one') do |tempfile|
      tempfile.puts '>one'
      tempfile.puts 'NNNNN'
      tempfile.close

      assert_equal "ID\tA\tC\n", `#{script_path} -w 5 -k 1 #{tempfile.path}`
    end
  end
  
  should 'give correct increments in window numbering' do
    Tempfile.open('one') do |tempfile|
      tempfile.puts '>one'
      tempfile.puts 'ATGCATGCAT' #10 letters long
      tempfile.close
      
      expected = "ID\tA\tC\n"+
      "one_0\t0.5\t0.5\n"+
      "one_1\t0.5\t0.5\n"+
      "one_leftover_2\t1.0\t0.0\n"

      assert_equal expected, `#{script_path} -w 4 -k 1 -m 2 #{tempfile.path}`
    end
  end
end
