#ifndef GENOME_H_
#define GENOME_H_

namespace structure { 

struct Chromosome { 
  typedef std::string orf_t;
  typedef std::pair<orf_t, int> gene_t;

  Chromosome () 
  : isCircular(false) 
  { 
  }

  Chromosome (const std::list<gene_t> & path, bool is_circ) 
  : isCircular(is_circ) 
  {
    size_t offset = 0; 
    for (const auto & gene : path) { 
      main_chromosome.insert(std::make_pair(offset++, gene));
    }  
  }

  inline void insert(const orf_t& gene, size_t offset, int sign) { 
    auto gen = gene_t(gene, sign);
    main_chromosome.insert(std::make_pair(offset, gen));
  } 

  inline void do_circular() { 
	isCircular = true; 
  } 

  inline std::map<size_t, gene_t>::const_iterator begin() const { 
    return main_chromosome.cbegin();	
  } 

  inline std::map<size_t, gene_t>::const_iterator end() const { 
    return main_chromosome.cend();
  } 

  inline std::map<size_t, gene_t>::const_reverse_iterator crbegin() const { 
    return main_chromosome.crbegin();	
  } 

  inline std::map<size_t, gene_t>::const_reverse_iterator crend() const { 
    return main_chromosome.crend();
  } 

  inline size_t size() const {
    return main_chromosome.size();      
  } 

  inline bool is_circular() const { 
    return isCircular;
  }
private: 
  bool isCircular; 
  std::map<size_t, gene_t> main_chromosome;
}; 

struct Genome {
  typedef std::string orf_t;
  typedef std::pair<std::string, size_t> coord_t; // (contig, offset)
  typedef std::pair<orf_t, int> gene_t;
  typedef structure::Chromosome chromosome_t; 
  
  Genome() {
  } 

  Genome(const std::string& nm) 
  : name(nm) 
  {
  } 

  inline void insert(const orf_t& gene, const std::string& chromosome, size_t offset, int sign, size_t start, size_t end) {
    const coord_t& p = std::make_pair(chromosome, offset);
    const gene_t& orf = std::make_pair(gene, sign);   
    main_genome[chromosome].insert(gene, offset, sign);
  }   

  inline void insert(const std::string& name_chr, const chromosome_t& chr) {
    main_genome.insert(std::make_pair(name_chr, chr)); 
  } 
	
  inline void registrate_circular_chr(std::string name) { 
    main_genome[name].do_circular();
  }  

  inline size_t size() const { 
    size_t sz = 0; 
    std::for_each(main_genome.begin(), main_genome.end(), [&](const std::pair<std::string, chromosome_t> chromosome) {
      sz += chromosome.second.size(); 
    });
    return sz; 
  }

  inline size_t count_chromosome() const { 
    return main_genome.size();
  } 

  inline std::string get_name() const {  
    return name;
  } 

  inline std::map<std::string, chromosome_t>::const_iterator begin() const { 
    return main_genome.cbegin();	
  } 

  inline std::map<std::string, chromosome_t>::const_iterator end() const { 
    return main_genome.cend();
  } 
private: 
  std::string name;
  std::map<std::string, chromosome_t> main_genome;
};

} 
#endif
