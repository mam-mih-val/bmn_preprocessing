//
// Created by mikhail on 4/20/22.
//

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

#include <boost/program_options.hpp>

#include <AnalysisTree/TaskManager.hpp>

#include "tracks_processor.h"

int main(int argv, char **argc){
  namespace po = boost::program_options;
  if( argv < 2 ){
    throw std::runtime_error( "No arguments are specified. Try ./run_qa --help for more information" );
  }
  std::vector<std::string> in_file_lists;
  std::vector<std::string> in_tree_names{ std::vector{std::string("rTree")} };
  std::string output_file{"output.root"};
  std::string output_tree{"eTree"};
  int n_events{-1};
  po::options_description options("Options");
  options.add_options()
      ("help,h", "Help screen")
      ("MC,m", "MC data sample")
      ("input,i", po::value<std::vector<std::string>>(&in_file_lists),"Input file list")
      ("tree-name,t", po::value<std::vector<std::string>>(&in_tree_names),"Input file list")
      ("output-tree-name", po::value<std::string>(&output_tree),"Output tree name")
      ("output,o", po::value<std::string>(&output_file),"Name of output file")
      ("events,N", po::value<int>(&n_events),"Number of analysing events");
  po::variables_map vm;
  po::parsed_options parsed = po::command_line_parser(argv, argc).options(options).run();
  po::store(parsed, vm);
  po::notify(vm);
  if (vm.count("help")){
    std::cout << options << std::endl;
    return 0;
  }

  auto* man = AnalysisTree::TaskManager::GetInstance();
  auto* task = new AnalysisTree::TracksProcessor();
  man->AddTask(task);

  man->Init(in_file_lists, in_tree_names);
  man->SetOutputName(output_file, output_tree);
  man->Run(n_events);
  man->Finish();

  return 0;
}