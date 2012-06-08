
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
      # Return a hash of Strings to 0, for each kmer of length k. For instance
      # empty_full_kmer_hash(1) => {'A'=>0, 'T'=>0, 'C'=>0, 'G'=>0}
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
      
      # Take a kmer hash, and merge those keys to the lowest lexigraphical form
      # (See Bio::Sequence::NA#lowest_lexigraphical_form for what this means)
      # When 2 keys are reverse complements they get merged into one hash entry, 
      # where the key is the lowest_lexigraphical_form of the two and the
      # value is the sum of the original 2 values
      # 
      # For instance {'A'=>2,'T'=>5} #=> {'A'=>7}  
      def self.merge_down_to_lowest_lexigraphical_form(hash)
        keys = empty_full_kmer_hash.keys
        
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