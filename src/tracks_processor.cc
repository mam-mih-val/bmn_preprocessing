//
// Created by mikhail on 11/30/20.
//

#include "tracks_processor.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TaskManager.hpp"

#include <AnalysisTree/DataHeader.hpp>
#include <random>
#include <TH2F.h>

namespace AnalysisTree {

void TracksProcessor::Init() {
  auto man = TaskManager::GetInstance();
  auto chain = man->GetChain();

  AddInputBranch("SimParticles");
  AddInputBranch("GlobalTracks");
  AddInputBranch("RecEventHeader");

//  AddInputBranch("SimParticles2GlobalTracks");
  in_sim_particles_ = chain->GetBranch("SimParticles");
  in_sim_particles_.Freeze();
  in_tracks_ = chain->GetBranch("GlobalTracks");
  in_tracks_.Freeze();
  in_event_header_ = chain->GetBranch("RecEventHeader");
  in_event_header_.Freeze();
  sim_particles_2_global_tracks_ = chain->GetMatching( "SimParticles", "GlobalTracks" );

  auto out_event_header_conf = BranchConfig("RecEventHeaderExt", DetType::kEventHeader);
  out_event_header_conf.AddField<float>("centrality", "(%) Based on multiplicity");
  out_event_header_ = Branch(out_event_header_conf);
  out_event_header_.SetMutable();

  auto in_sim_particles_conf = chain->GetConfiguration()->GetBranchConfig("SimParticles");
  auto out_sim_particles_conf = in_sim_particles_conf.Clone("SimParticlesExt", DetType::kParticle);
  out_sim_particles_conf.AddField<float>("Ekin", "kinetic energy");
  out_sim_particles_conf.AddField<float>("y_cm", "center-of-mass rapidity");

  out_sim_particles_ = Branch(out_sim_particles_conf);
  out_sim_particles_.SetMutable();

  auto in_tracks_conf = chain->GetConfiguration()->GetBranchConfig("GlobalTracks");
  auto out_tracks_conf = in_tracks_conf.Clone("GlobalTracksExt", DetType::kParticle);
  out_tracks_conf.AddField<bool>("has_tof_hit", "has either TOF400 or TOF700 hit");
  out_tracks_conf.AddField<bool>("is_primary", "is primary particle");
  out_tracks_conf.AddField<float>("y_cm", "center-of-mass rapidity");
  out_tracks_conf.AddField<float>("efficiency", "efficiency");
  out_tracks_conf.AddField<float>("tof_efficiency", "efficiency in TOF acceptance");
  out_tracks_conf.AddField<float>("weight", "efficiency > 0.01 ? 1/efficiency : 0.0");
  out_tracks_conf.AddField<float>("tof_weight", "tof_efficiency > 0.01 ? 1/efficiency : 0.0");
  out_tracks_conf.AddField<float>("x_fhcal", "Linear extrapolation of track coordinates to FHCal plane");
  out_tracks_conf.AddField<float>("y_fhcal", "Linear extrapolation of track coordinates to FHCal plane");
  out_tracks_conf.AddField<bool>("in_fhcal", "If the particle is in FHCal acceptance");

  out_tracks_ = Branch(out_tracks_conf);
  out_tracks_.SetMutable();

  man->AddBranch(&out_sim_particles_);
  man->AddBranch(&out_tracks_);
  man->AddBranch(&out_event_header_);

  FHCalQA();
  ReadEfficiency();
}

void TracksProcessor::Exec() {
  using AnalysisTree::Particle;
  this->FillEventHeader();
  this->LoopRecTracks();
  if (is_mc_) {
    this->LoopSimParticles();
  }
}

void TracksProcessor::LoopRecTracks() {
  using AnalysisTree::Particle;
  out_tracks_.ClearChannels();
  auto field_out_pT = out_tracks_.GetField("pT");
  auto field_out_ycm = out_tracks_.GetField("y_cm");
  auto field_out_efficiency = out_tracks_.GetField("efficiency");
  auto field_out_weight = out_tracks_.GetField("weight");
  auto field_out_tof_efficiency = out_tracks_.GetField("tof_efficiency");
  auto field_out_tof_weight = out_tracks_.GetField("tof_weight");
  auto field_out_beta400 = out_tracks_.GetField("beta400");
  auto field_out_beta700 = out_tracks_.GetField("beta700");
  auto field_out_has_tof_hit = out_tracks_.GetField("has_tof_hit");
  auto field_out_is_primary = out_tracks_.GetField("is_primary");
  auto field_out_pid = out_tracks_.GetField("pid");
  auto field_out_mass = out_tracks_.GetField("mass");
  auto field_out_rapidity = out_tracks_.GetField("rapidity");

  auto field_x_fhcal = out_tracks_.GetField("x_fhcal");
  auto field_y_fhcal = out_tracks_.GetField("y_fhcal");
  auto field_in_fhcal = out_tracks_.GetField("in_fhcal");

  auto field_x_last = out_tracks_.GetField("x_last");
  auto field_y_last = out_tracks_.GetField("y_last");
  auto field_z_last = out_tracks_.GetField("z_last");
  auto field_tx_last = out_tracks_.GetField("tx_last");
  auto field_ty_last = out_tracks_.GetField("ty_last");

  auto field_sim_pid = in_sim_particles_.GetField("pid");
  auto field_sim_mass = in_sim_particles_.GetField("mass");
  auto field_sim_mother_id = in_sim_particles_.GetField("mother_id");

  auto y_beam = data_header_->GetBeamRapidity();

  for (size_t i=0; i<in_tracks_.size(); ++i) {
    auto in_track = in_tracks_[i];
    auto sim_particle_idx = sim_particles_2_global_tracks_->GetMatchInverted( i );
    if( sim_particle_idx < 0 )
      continue;
    auto sim_particle = in_sim_particles_[sim_particle_idx];
    auto pid = (int) sim_particle[field_sim_pid];
    auto mother_id = (int) sim_particle[field_sim_mother_id];
    auto mass = sim_particle[field_sim_mass];
    auto out_particle = out_tracks_.NewChannel();
    out_particle.CopyContent( in_track );
    out_particle.SetValue( field_out_mass, (float) mass );
    out_particle.SetValue( field_out_pid, (int) pid );
    out_particle.SetValue( field_out_is_primary, mother_id == -1 );
    auto rapidity = out_particle[field_out_rapidity];
    out_particle.SetValue( field_out_ycm, (float)rapidity - y_beam );
    auto beta400 = out_particle[field_out_beta400];
    auto beta700 = out_particle[field_out_beta700];
    bool has_tof_hit = beta400 > -990. || beta700 > -990.;
    out_particle.SetValue( field_out_has_tof_hit, has_tof_hit );
    auto pT = out_particle[field_out_pT];
    auto ycm = rapidity - y_beam;
    auto [efficiency, tof_efficiency] = FindEfficiency( pid, pT, ycm );
    double weight = efficiency > 0.01 ? 1.0/efficiency : 0.0;
    double tof_weight = tof_efficiency > 0.01 ? 1.0/tof_efficiency : 0.0;

    out_particle.SetValue( field_out_efficiency, float(efficiency) );
    out_particle.SetValue( field_out_weight, float(weight) );

    out_particle.SetValue( field_out_tof_efficiency, float(tof_efficiency) );
    out_particle.SetValue( field_out_tof_weight, float(tof_weight) );

    auto x = out_particle[field_x_last];
    auto y = out_particle[field_y_last];
    auto z = out_particle[field_z_last];
    auto tx = out_particle[field_tx_last];
    auto ty = out_particle[field_ty_last];
    auto [x_fhcal, y_fhcal] = ProjectToFHCalPlane( x, y, z, tx, ty );
    out_particle.SetValue(field_x_fhcal, float(x_fhcal));
    out_particle.SetValue(field_y_fhcal, float(y_fhcal));
    auto in_fhcal = true;
    if( x_fhcal < -30 )
      in_fhcal = false;
    if( x_fhcal > 140 )
      in_fhcal = false;
    if( fabs(y_fhcal) > 60  )
      in_fhcal = false;
    out_particle.SetValue(field_in_fhcal, bool(in_fhcal));
  }
}
void TracksProcessor::LoopSimParticles() {
  using AnalysisTree::Particle;
  out_sim_particles_.ClearChannels();
  auto field_out_ycm = out_sim_particles_.GetField("y_cm");
  auto field_out_Ekin = out_sim_particles_.GetField("Ekin");
  auto field_out_mass = out_sim_particles_.GetField("mass");
  auto field_out_p = out_sim_particles_.GetField("p");
  auto field_out_rapidity = out_sim_particles_.GetField("rapidity");

  auto y_beam = data_header_->GetBeamRapidity();

  for (size_t i=0; i<in_sim_particles_.size(); ++i) {
    auto in_particle = in_sim_particles_[i];
    auto out_particle = out_sim_particles_.NewChannel();
    out_particle.CopyContent(in_particle);
    auto rapidity = out_particle[field_out_rapidity];
    out_particle.SetValue( field_out_ycm, (float)rapidity - y_beam );
    auto p = out_particle[field_out_p];
    auto mass = out_particle[field_out_mass];
    auto Ekin = sqrt(p*p + mass*mass) - mass;
    out_particle.SetValue(field_out_Ekin, (float) Ekin);
  }
}
void TracksProcessor::FHCalQA() {
  auto man = TaskManager::GetInstance();
  auto data_header = man->GetDataHeader();
  auto module_pos = data_header->GetModulePositions( 0 );
  auto h2_module_positions_ = new TH2F("fhcal_module_positions", ";x (cm);y (cm)", 500, -100, 100, 500, -100, 100);
  for( int idx = 0; idx < 54; idx++ ){
    auto module = module_pos.Channel(idx);
    h2_module_positions_->Fill(module.GetX(), module.GetY(), module.GetId() );
  }
  h2_module_positions_->Write();
  auto h2_left_module_positions_ = new TH2F("scwall_module_positions", ";x (cm);y (cm)", 500, -150, 150, 500, -100, 100);
  for( int idx = 54; idx < 54+174; idx++ ){
    auto module = module_pos.Channel(idx);
    h2_left_module_positions_->Fill(module.GetX(), module.GetY(), module.GetId() );
  }
  h2_left_module_positions_->Write();
}

void TracksProcessor::FillEventHeader() {
  auto field_in_M = in_event_header_.GetField("M");
  auto field_out_centrality = out_event_header_.GetField("centrality");
  auto multiplicity = in_event_header_.GetDataRaw<EventHeader*>()->GetField<int>(field_in_M.GetFieldId());
  auto centrality = -1.0f;
  int idx = 0;
  float bin_edge = multiplicity_edges_[idx];
  while( multiplicity < bin_edge &&
         idx < multiplicity_edges_.size()-1 ){
    idx++;
    bin_edge = multiplicity_edges_[idx];
  }
  centrality = (centrality_percentage_[idx-1] + centrality_percentage_[idx])/2.0f;
  out_event_header_.GetDataRaw<EventHeader*>()->SetField(centrality, field_out_centrality.GetFieldId());
}
void TracksProcessor::ReadEfficiency() {
  if( efficiency_file_name_.empty() )
    return;
  efficiency_file_ = TFile::Open(efficiency_file_name_.c_str());
  if( !efficiency_file_ ) {
    std::cerr << "No such file: " << efficiency_file_name_ << std::endl;
    return;
  }
  // Efficiency of tracking
  efficiency_file_->GetObject("efficiency_2212", efficiency_2212_);
  if( !efficiency_2212_ )
    std::cerr << "File " << efficiency_file_name_ << " does not contain histogram for proton efficiency" << std::endl;
  efficiency_file_->GetObject("efficiency_-211", efficiency_m211_);
  if( !efficiency_m211_ )
    std::cerr << "File " << efficiency_file_name_ << " does not contain histogram for negative pion efficiency" << std::endl;
  // Efficiency of tracking in TOF-acceptance
  efficiency_file_->GetObject("efficiency_2212_tof", efficiency_2212_tof_);
  if( !efficiency_2212_tof_ )
    std::cerr << "File " << efficiency_file_name_ << " does not contain histogram for proton TOF efficiency" << std::endl;
  efficiency_file_->GetObject("efficiency_-211_tof", efficiency_m211_tof_);
  if( !efficiency_m211_tof_ )
    std::cerr << "File " << efficiency_file_name_ << " does not contain histogram for negative pion TOF efficiency" << std::endl;
}
std::tuple<double, double> TracksProcessor::FindEfficiency(int pid, double pT, double y) {
  TH2F* hist{nullptr};
  TH2F* hist_tof{nullptr};
  if( pid == 2212 ) {
    hist = efficiency_2212_;
    hist_tof = efficiency_2212_tof_;
  }
  if( pid == -211 ) {
    hist = efficiency_m211_;
    hist_tof = efficiency_m211_tof_;
  }
  double eff = 1.0;
  double eff_tof = 1.0;
  if( hist ) {
    auto y_bin = hist->GetXaxis()->FindBin(y);
    auto pT_bin = hist->GetYaxis()->FindBin(pT);
    eff = hist->GetBinContent(y_bin, pT_bin);
  }
  if( hist_tof ) {
    auto y_bin = hist_tof->GetXaxis()->FindBin(y);
    auto pT_bin = hist_tof->GetYaxis()->FindBin(pT);
    eff_tof = hist_tof->GetBinContent(y_bin, pT_bin);
  }
  return {eff, eff_tof};
}
std::tuple<double, double>
TracksProcessor::ProjectToFHCalPlane(double x, double y, double z,
                                     double Tx, double Ty) {
  auto dz = 900. - z;
  auto dx = Tx*dz;
  auto dy = Ty*dz;
  auto x_fhcal = x + dx;
  auto y_fhcal = y + dy;
  return {x_fhcal, y_fhcal};
}
} // namespace AnalysisTree