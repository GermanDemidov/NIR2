#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include "mbgraph_history.h"
#include "genome_match.h" //FIXME REMOVE LATER
#include "writer/Wstats.h"
#include "writer/Wdots.h"

template<class graph_t>
struct Algorithm { 
  Algorithm(const std::shared_ptr<graph_t>& gr, size_t size_component, size_t rds) 
  : graph(gr) 
  , canformQoo(true)
  , rounds(rds)
  , max_size_component(size_component)
  , write_stats("stats.txt") 
  {
    //generate_disjoint_colors();
  } 

  void convert_to_identity_bgraph(const ProblemInstance<structure::Mcolor>& cfg);

  edges_t get_bad_edges() const; 

  //Stage 1: process insertion/deletion events. Balanced graph
  bool stage3();

private: 
  typedef structure::Mcolor mcolor_t; 
  typedef structure::Mularcs<mcolor_t> mularcs_t; 
  typedef event::TwoBreak<mcolor_t> twobreak_t;
  typedef event::FakeTwoBreak<mcolor_t> fake_twobreak_t;        
  typedef event::InsDel<mcolor_t> insertion_t;
  typedef event::TandemDuplication<mcolor_t> tandem_duplication_t;
	
  //Stage 2: Simple paths  
  bool stage1(size_t number_splits); 
  size_t process_simple_path(path_t& path, size_t number_splits);	

  //Stage 2: process non-mobile edges
  bool stage2();
  bool stage22();
  bool is_mobility_edge(const vertex_t& x, const vertex_t& y) const;
  bool is_mobility_edge(const vertex_t& x, const mcolor_t& color, const vertex_t& y) const;
  size_t is_mobility_edge_score(const vertex_t& x, const mcolor_t& color, const vertex_t& y) const;

  //Stage -: process tandem duplication events
  bool stage4_td(); 
  bool stage4_rtd(); 

  //Stage 3: process components, 
  //utility::equivalence<vertex_t> split_on_components(std::map<vertex_t, std::map<mcolor_t, std::set<arc_t> > >& EC, std::map<vertex_t, std::map<mcolor_t, std::set<arc_t> > >& EI, const std::set<mcolor_t>& colors);
  utility::equivalence<vertex_t> split_on_components(std::map<vertex_t, std::set<arc_t> >& EC, std::map<vertex_t, std::set<arc_t> >& EI, const mcolor_t& Q);
  bool stage5_1(); 
  bool stage5_2();

  //Stage 4: Clone approach
  bool stage7();
        
  //Stage 10: convert duplication to tandem duplication   
  bool stage4_conv_to_td();

  //Stage 99: brutoforce stage
  bool stage6();
  size_t calculate_cost(const vertex_t& y, const mularcs_t& mularcs_x, const mularcs_t& mulacrs_y); 
  std::set<arc_t> create_minimal_matching(const std::set<vertex_t>& vertex_set); 
  size_t process_minimal_matching(const std::set<arc_t>& matchings);

  size_t check_postponed_deletions() const;
  void generate_disjoint_colors();
private: 
  std::shared_ptr<graph_t> graph; 

  bool canformQoo;  // safe choice, at later stages may change to false
  const size_t rounds;
  const size_t max_size_component; 

  edges_t insertions;
  edges_t postponed_deletions; 
  std::map<std::pair<vertex_t, mcolor_t>, vertex_t> mother_verteces;
 
  std::set<edge_t> clone_edges; 
  std::set<std::set<mcolor_t> > disjoint_subset_colors;  

  writer::Wstats write_stats;
  writer::Wdots<graph_t, ProblemInstance<mcolor_t> > write_dots; 
};

template<class graph_t>
void Algorithm<graph_t>::convert_to_identity_bgraph(const ProblemInstance<structure::Mcolor>& cfg) {
  std::unordered_set<size_t> print_dots;
  size_t stage = 0; 
  bool isChanged = false;
  size_t count_changes = 0;
  bool process_compl = true; 
  auto saveInfoLambda = [&](size_t st) -> void { 
    if ((print_dots.count(st) == 0) && !isChanged) {
      print_dots.insert(st);
      Statistics<graph_t> stat(graph); 
      write_stats.print_all_statistics(st, stat, *graph);
      write_dots.save_dot(*graph, cfg, st);
    } 
  };
  saveInfoLambda(0);
  size_t brf = 0;
 
#ifdef VERSION2
  if (cfg.get_stages() >= 1) { 
    std::cerr << "Stage: 1 (indel stage)" << std::endl;
    stage3();
    saveInfoLambda(1);
  }
#endif

  isChanged = true;
  while(isChanged) {
    isChanged = false; 
    stage = 2; 

#ifdef VERSION2
    
    for (size_t i = 1; i <= rounds && !isChanged; ++i) {   
      std::cerr << "Rounds " << i << std::endl;

      graph->update_number_of_splits(i);
    
      if ((cfg.get_stages() >= 2) && !isChanged) {
        std::cerr << "Stage: 2 Good path " << stage << std::endl;

        isChanged = stage1(i);	
 
        if (!isChanged) {
          std::cerr << "Stage: 2 Non-mobile edges " << stage << std::endl;
	  isChanged = stage22();
        }

        saveInfoLambda(stage++);
      }

      if ((cfg.get_stages() >= 3) && !isChanged) { // STAGE 4, somewhat unreliable
        std::cerr << "Stage: 3 Split on components " << stage << std::endl;
        isChanged = stage5_1(); // cut the graph into connected components
        saveInfoLambda(stage++); 
      }

      if ((cfg.get_stages() >= 4) && !isChanged) {
        std::cerr << "Stage: 4 Clone approach " << stage << std::endl;
        isChanged = stage7();
        saveInfoLambda(stage++);
      }
    } 

    if (isChanged) {
      count_changes = 0;
    }
 
    if (!isChanged && count_changes != 2) { 
      std::cerr << "Change canformQoo " << count_changes << std::endl;
      canformQoo = !canformQoo; // more flexible
      ++count_changes;
      isChanged = true;
    }

    /*if ((max_size_component != 0) && !isChanged) {
      std::cerr << "Brute force stage: " << std::endl;
      graph->update_number_of_splits(3);
      isChanged = stage6();
      saveInfoLambda(stage++);
    }*/

#else
   stage = 2; 

   graph->update_number_of_splits(1);
    
    if ((cfg.get_stages() >= 1) && !isChanged) {
      //std::cerr << "Stage: 1" << std::endl;
      isChanged = stage1(1);	
      saveInfoLambda(stage++);
    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
      //std::cerr << "Stage: 2" << std::endl;
      isChanged = stage2();
      saveInfoLambda(stage++);
    }

    if ((cfg.get_stages() >= 3) && !isChanged) { // STAGE 3, somewhat unreliable
      //std::cerr << "Stage: 3" << std::endl;
      isChanged = stage5_1(); // cut the graph into connected components
      
      if (!isChanged) { 
         isChanged = stage5_2(); // process 4-cycles
      }  
     
      if (canformQoo && !isChanged) {
	isChanged = true;
	canformQoo = false; // more flexible
      }    

      saveInfoLambda(stage++);
    }

    if ((cfg.get_stages() >= 4) && !isChanged) {
      //std::cerr << "Stage: 4" << std::endl;
      graph->update_number_of_splits(graph->count_local_graphs());
      isChanged = stage2();
      saveInfoLambda(stage++);
    }

    if (process_compl && !cfg.get_completion().empty() && !isChanged) {     
      //std::cerr << "Manual Completion Stage" << std::endl;
      auto completion = cfg.get_completion();
      for(auto il = completion.begin(); il != completion.end(); ++il) {
	graph->apply_two_break(*il);
      }

      process_compl = false;
      isChanged = true;
    }
#endif
  }	

  write_dots.save_dot(*graph, cfg, 99);

  graph->change_history();

  size_t bad_postponed_deletions = check_postponed_deletions();
  if (bad_postponed_deletions != 0) {
    std::cerr << "We have problem with " << bad_postponed_deletions << " edges, corresponding postponed deletions." << std::endl;
    std::cerr << "If you have indentity breakpoint graph after stages, please contact us." << std::endl;
    exit(1);
  } 

  write_stats.print_history_statistics(*graph, get_bad_edges());
}          

template<class graph_t>
void Algorithm<graph_t>::generate_disjoint_colors() {
  for (auto it = graph->cbegin_T_consistent_color(); it != graph->cend_T_consistent_color(); ++it) {
    disjoint_subset_colors.insert(std::set<mcolor_t>({*it}));
    /*for (auto ic = graph->cbegin_T_consistent_color(); ic != graph->cend_T_consistent_color(); ++ic) {
      mcolor_t inter_color(*it, *ic, mcolor_t::Intersection);
      if (inter_color.empty()) {
        disjoint_subset_colors.insert(std::set<mcolor_t>({*it, *ic}));
      } 
    }*/
  } 
} 

template<class graph_t>
edges_t Algorithm<graph_t>::get_bad_edges() const { 
  edges_t answer; 
  for (const auto &edge: insertions) { 
    answer.insert(edge.first, edge.second);
  } 

  for (const auto &edge: postponed_deletions) {
    answer.insert(edge.first, edge.second);
  }
  return answer;
} 

#include "Stage1.h" 
#include "Stage2.h"
#include "Stage3.h"
//#include "Stage4.h"
#include "Stage5.h"
#include "Stage6.h"
#include "Stage7.h"

#endif
