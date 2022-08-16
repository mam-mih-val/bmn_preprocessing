//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

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
  TracksProcessor();
  ~TracksProcessor();
  void Init() override;
  void Exec() override;
  void Finish() override {}
  void SetEfficiencyFileName(const std::string &efficiency_file_name) {
    efficiency_file_name_ = efficiency_file_name;
  }

protected:
  void InitFields();
  void LoopRecTracks();
  void LoopSimParticles();
  void LoopScWallModules();
  void FillEventHeader();
  void FHCalQA();
  void ReadEfficiency();
  std::tuple<double, double> FindEfficiency(int pid, double pT, double y);
  std::tuple<double, double> ProjectToFHCalPlane(double x, double y, double z, double Tx, double Ty);

private:
  bool is_mc_ = true;

  Branch in_tracks_;
  Branch in_sim_particles_;
  Branch in_event_header_;
  Branch in_fhcal_modules_;

  Branch out_tracks_;
  Branch out_sim_particles_;
  Branch out_event_header_;
  Branch out_scwall_hits_;

  std::vector<float> centrality_percentage_{0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100 };
  std::vector<int> multiplicity_edges_{ 249, 155, 129, 108, 90, 74, 60, 49, 39, 24, 14, 7, 2, 1, 0 };

  std::string efficiency_file_name_;
  TFile* efficiency_file_;
  TH2F* efficiency_2212_;
  TH2F* efficiency_2212_tof_;
  TH2F* efficiency_m211_;
  TH2F* efficiency_m211_tof_;

  Matching* sim_particles_2_global_tracks_;
  ModulePositions module_positions_;

  struct scwall_fields;
  struct fhcal_fields;
  std::unique_ptr<scwall_fields> scwall_fields_;
  std::unique_ptr<fhcal_fields> fhcal_fields_;
};

}
#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
