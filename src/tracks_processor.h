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

private:
  bool is_mc_ = true;

  Branch in_tracks_;
  Branch in_sim_particles_;

  Branch out_tracks_;
  Branch out_sim_particles_;

  Matching* sim_particles_2_global_tracks_;
};

}
#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
