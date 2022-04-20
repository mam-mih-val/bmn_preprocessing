#add_executable(rapidity tracks_processor.cc tracks_processor.h sim_particles_processor.cc sim_particles_processor.h)
add_executable(pre_process tracks_processor.cc tracks_processor.h main.cc)
target_link_libraries(pre_process ROOT::EG AnalysisTreeBase AnalysisTreeInfra Boost::program_options)