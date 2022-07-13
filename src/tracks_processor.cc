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

  out_tracks_ = Branch(out_tracks_conf);
  out_tracks_.SetMutable();

  man->AddBranch(&out_sim_particles_);
  man->AddBranch(&out_tracks_);
  man->AddBranch(&out_event_header_);

  FHCalQA();
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
  auto field_out_ycm = out_tracks_.GetField("y_cm");
  auto field_out_beta400 = out_tracks_.GetField("beta400");
  auto field_out_beta700 = out_tracks_.GetField("beta700");
  auto field_out_has_tof_hit = out_tracks_.GetField("has_tof_hit");
  auto field_out_is_primary = out_tracks_.GetField("is_primary");
  auto field_out_pid = out_tracks_.GetField("pid");
  auto field_out_mass = out_tracks_.GetField("mass");
  auto field_out_rapidity = out_tracks_.GetField("rapidity");

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
} // namespace AnalysisTree