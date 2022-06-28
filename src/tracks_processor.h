//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <TFile.h>
#include <TTree.h>

#include <AnalysisTree/AnalysisTask.hpp>

#include <AnalysisTree/Matching.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Particle.hpp>
#include <AnalysisTree/BranchConfig.hpp>
#include <AnalysisTree/EventHeader.hpp>

#include <memory>
#include <string>

namespace AnalysisTree {

class TracksProcessor : public Task {

public:
  void Init() override;
  void Exec() override;
  void Finish() override {}

protected:
  void LoopRecTracks();
  void LoopSimParticles();
  void FillEventHeader();
  void FHCalQA();

private:
  bool is_mc_ = true;

  Branch in_tracks_;
  Branch in_sim_particles_;
  Branch in_event_header_;

  Branch out_tracks_;
  Branch out_sim_particles_;
  Branch out_event_header_;

  std::vector<float> centrality_percentage_{0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100 };
  std::vector<int> multiplicity_edges_{ 249, 155, 129, 108, 90, 74, 60, 49, 39, 24, 14, 7, 2, 1, 0 };

  Matching* sim_particles_2_global_tracks_;
};

}
#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
