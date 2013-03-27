require 'helper'
require 'tempfile'
require 'open3'

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

  should "print help when no arguments are given" do
    command = "#{script_path}"
    Open3.popen3(command) do |stdin, stdout, stderr|
      assert stderr.readlines[0].match(/^Usage: kmer_counter/)
    end
  end

  should 'work with lowercase' do
    Tempfile.open('one') do |tempfile|
      tempfile.puts '>one'
      tempfile.puts 'acagt'
      tempfile.close

      assert_equal "ID\tA\tC\none_0\t0.6\t0.4\n", `#{script_path} -w 5 -k 1 #{tempfile.path}`
    end
  end

  should 'by default count contigs greater than 2kb but less than 5kb' do
    Tempfile.open('one') do |tempfile|
      tempfile.puts '>one'
      tempfile.puts 'A'*2500
      tempfile.close

      assert_equal "ID\tA\tC\none_leftover_0\t1.0\t0.0\n", `#{script_path} -k 1 #{tempfile.path}`
    end
  end

  should 'by default count contigs greater than 2kb but less than 5kb' do
    Tempfile.open('one') do |tempfile|
      tempfile.puts '>one'
      tempfile.puts 'A'*7500
      tempfile.close

      assert_equal "ID\tA\tC\none_0\t1.0\t0.0\none_leftover_1\t1.0\t0.0\n", `#{script_path} -k 1 #{tempfile.path}`
    end
  end

  should 'work simulated example with kmer length = 2' do
    expected = %w(ID	AA	AC	AG	AT	CA	CC	CG	GA	GC	TA).join("\t")+"\n"+
    %w(random_leftover_0	0.1111111111111111	0.13131313131313133	0.1414141414141414	0.0707070707070707	0.1717171717171717	0.1111111111111111	0.020202020202020204	0.1414141414141414	0.050505050505050504	0.050505050505050504).join("\t")+"\n"

    assert_equal expected, `#{script_path} -k 2 -m 1 #{File.join(TEST_DATA_DIR,'100random.fa')}`
  end
end
