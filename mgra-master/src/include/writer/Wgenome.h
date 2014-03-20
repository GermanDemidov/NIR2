#ifndef WGENOME_H_
#define WGENOME_H_

namespace writer {

template <class genome_t>
struct Wgenome {
  void save_genomes(const std::vector<genome_t>& genomes, bool isEmptyTarget) const { 
    for (const auto &genome: genomes) {
      save_genome_in_text_format(genome, isEmptyTarget);
    } 
  } 

private: 
  void save_genome_in_text_format(const genome_t& genome, bool isEmptyTarget) const;
};

} 

template <class genome_t>
void writer::Wgenome<genome_t>::save_genome_in_text_format(const genome_t& genome, bool isEmptyTarget) const { 
  const std::string& outname = genome.get_name();  
  std::ofstream out((outname + ".gen").c_str());

  out << "# Genome " << outname << std::endl;

  std::string chr_title;

  if (isEmptyTarget) { 
    chr_title = "chromosome"; 
  } else { 
    chr_title = "CAR";
  } 

  size_t number_circular = 0; 
  size_t length_circular = 0;

  for(const auto &chromosome : genome) {
    out << std::endl;

    if (chromosome.second.is_circular()) {
      ++number_circular;
      length_circular += chromosome.second.size();
      out << "# circular ";
    } else {
      out << "# linear ";
    }

    out << chr_title << " of length " << chromosome.second.size() << " follows:" << std::endl;

    if (chromosome.second.begin()->second.second == 1 || (--chromosome.second.end())->second.second == 1) {
      for(const auto &block : chromosome.second) {
	out << ((block.second.second == -1)?"-":"+") << block.second.first << " ";
      }
    } else {
      for(auto ip = chromosome.second.crbegin(); ip != chromosome.second.crend(); ++ip) {
	out << ((ip->second.second == -1)?"+":"-") << ip->second.first << " ";
      }
    }

    if (chromosome.second.is_circular()) {			
      out << "@" << std::endl;
    } else { 
      out << "$" << std::endl;
    } 
  }

  out << std::endl << "# Reconstructed genome " << outname << " has " << genome.count_chromosome() << " " << chr_title << "s" << std::endl;
#ifdef VERSION2
  std::cout << std::endl << "# Reconstructed genome " << outname << " has " << genome.count_chromosome() << " " << chr_title << "s" << std::endl;
#endif

  if (number_circular) {
    out << "#\tout of which " << number_circular << " are circular of total length " << length_circular << std::endl;
#ifdef VERSION2
    std::cout << "#\tout of which " << number_circular << " are circular of total length " << length_circular << std::endl;
#endif
  }
  out.close();
}

#endif
