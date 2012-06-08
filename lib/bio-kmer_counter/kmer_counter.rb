# Initialise the hash of the different
module Bio
  class Sequence
    class NA
      # Return the current object or its reverse complement, whichever
      # has the sequence that comes first in lexigraphical (alphabetical)
      # order
      def lowest_lexigraphical_form
        rev = self.reverse_complement
        to_s < rev.to_s ? self : rev
      end
    end

    class Kmer
      def self.empty_full_kmer_hash(k=4)
        return @empty_full_hash.dup unless @empty_full_hash.nil?
        
        counts = {}
                
        ordered_possibilities = %w(A T C G)
        keys = ordered_possibilities
        (k-1).times do
          keys = keys.collect{|k| ordered_possibilities.collect{|n| "#{k}#{n}"}.flatten}.flatten
        end
        
        keys.each do |key|
          counts[key] = 0
        end
        counts
      end
      
      def self.merge_down_to_lowest_lexigraphical_form(hash)
        keys = empty_full_kmer_hash.keys
        
        # remove keys (kmers) that are not in their lowest lexigraphical form.
        # Because all possible kmers are already present, the reverse complement of these must already be present in the array
        keys.select! do |key|
          Bio::Sequence::NA.new(key).lowest_lexigraphical_form.to_s.upcase == key
        end
        
        new_hash = {}
        hash.each do |kmer, count|
          key = Bio::Sequence::NA.new(kmer).lowest_lexigraphical_form.to_s.upcase
          new_hash[key] ||= 0
          new_hash[key] += count
        end
        return new_hash
      end
    end
  end
end