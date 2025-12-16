#ifndef ANA_H__
#define ANA_H__

// Utility
#include <utility>
#include <vector>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH2.h>
#include <TGraph2D.h>
#include <TF2.h>
#include <cassert>
#include <sstream>
#include <string>
#include <TLorentzVector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

// Fun4All
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

// Event
#include <Event/Event.h>
#include <Event/packet.h>
#include <ffaobjects/EventHeaderv1.h>

//Trigger
#include <calotrigger/TriggerRunInfov1.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <calotrigger/LL1Out.h>
#include <calotrigger/LL1Outv1.h>
#include <calotrigger/TriggerPrimitive.h>
#include <calotrigger/TriggerPrimitivev1.h>
#include <calotrigger/TriggerPrimitiveContainer.h>
#include <calotrigger/TriggerPrimitiveContainerv1.h>
#include <calotrigger/TriggerDefs.h>

// Global vertex
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>

// G4
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

// Tower includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfoContainerSimv1.h>
#include <calobase/TowerInfoSimv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfov2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfov3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfov4.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/PhotonClusterv1.h>

// MBD
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtSimContainerV1.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdGeom.h>

// phool
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <zdcinfo/Zdcinfov2.h>
#include <epd/EpdGeom.h>



#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

// Centrality MB
#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <centrality/CentralityInfov1.h>

#include <ffarawobjects/Gl1Packet.h>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TH2F;
class TH1;
class Gl1Packet;
class PHG4Shower;
class TriggerPrimitive;
class TriggerPrimitiveContainer;
class LL1Out;
class TruthNeutralMeson;

class Ana : public SubsysReco
{
 public:
  //! constructor
  Ana(const std::string &name = "Ana",  const char* outname="DST-00021615-0000.root", bool _isMC=false);

  //! destructor
  virtual ~Ana();

  //! full initialization
  int Init(PHCompositeNode *);
  void InitOutputFile();
  void InitTree();
  
  //! event processing method
  int process_event(PHCompositeNode *);
  int ProcessGlobalEventInfo(PHCompositeNode *);
  void ProcessFillTruthMesonParticle(PHCompositeNode *, float truthvz);
  void ProcessFillTruthPhotonParticle(float truthvz);

  //! end of run method
  int End(PHCompositeNode *);

  Double_t GetShiftedEta(float _vz, float _eta);

  int process_g4hits(PHCompositeNode *);
  int process_g4cells(PHCompositeNode *);
  int process_towers(PHCompositeNode *);
  int process_clusters(PHCompositeNode *);
  bool FindConversion(PHG4TruthInfoContainer *, int trackid, float energy);
  int Getpeaktime(TH1 *h);
  PHG4Particle* g4_to_top(PHG4Particle* p);
  bool is_hadron(int pdg){return (std::abs(pdg) >= 23); }
  HepMC::GenParticle* get_hepmc_particle(int barcode);
  int find_hepmc_pdg(HepMC::GenParticle* p);
  bool is_g4_hadron(PHG4Particle* p);


  void SetMbdZVtxCut(float _zvtxcut)
  {
    if(_zvtxcut<=0) IsVtxCut = false;
    else IsVtxCut = true;
    vzcut = _zvtxcut;
  };
  
  void SetRunNumber(int RunNumber)
  {
    m_run_number = RunNumber;
  };

  void SetScaledowns(int scaledowns[])
  {
    for(int i=0; i<64;i++){
      m_scaledowns[i] = scaledowns[i];
      std::cout << "prescale set " << i << " : " << m_scaledowns[i] << std::endl;
    }
  };

  void SetTriggerEmulator(bool _useEmulator)
  {
    useEmulator = _useEmulator;
  };

  void Detector(const std::string &name) { detector = name; }
  void ProcessClearBranchVar();

 protected:
  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm = nullptr;
  TFile *outfile = nullptr;
  TTree *towerntuple = nullptr;
  TH1D * profilehist;
  TH1D * histlumicount;
  Long64_t mbdlivecount=0;
  std::vector<GlobalVertex::VTXTYPE> m_vertex_type_sel{GlobalVertex::UNDEFINED};
  std::unordered_map<int,double> mbd_t0_corr;

  float m_vertex;
  int m_run_number=-999;
  short m_scaledowns[64];
  float _mbd_charge_threshold = 0.4;
  short m_mbd_hits_south = 0;
  short m_mbd_hits_north = 0;
  float m_mbd_mean_time = -999;
  float m_mbd_time_north = -999;
  float m_mbd_time_south = -999;
  float m_mbd_sumcharge_south = 0;
  float m_mbd_sumcharge_north = 0;

  float m_zdc_energy_south[4];
  float m_zdc_energy_north[4];
  float m_zdcinfo_sum_energy_south;
  float m_zdcinfo_sum_energy_north;
  float m_zdcinfo_vtx;
  std::string m_zdc_node_name= "TOWERINFO_CALIB_SEPD"; 
  std::string m_zdcinfo_node_name= "Zdcinfo"; 

  float m_sepd_energy[768];
  int m_sepd_arm[768];
  float m_sepd_phi[768];
  float m_sepd_radius[768];
  std::string m_sepd_node_name = "TOWERINFO_CALIB_SEPD"; 

  float vx, vy, vz;
  float vzcut;
  bool IsVtxCut = true;
  bool ScaledTriggerBit[64];
  bool LiveTriggerBit[64];

  int count=0;

  bool isMC;
  int processId;
  static const int nParticleTruth = 10000;
  float truth_vz=-999;
  float truth_vx=-999;
  float truth_vy=-999;

  bool doEMCal = false;
  bool doHCal = false;
  bool doMBD = false;
  bool doZDC = false;
  bool isPythia = false;
  bool isSingleGun = false;

  TriggerAnalyzer *triggeranalyzer{nullptr};
  bool useEmulator = false;

  PHHepMCGenEventMap* m_genevtmap = nullptr;
  PHG4TruthInfoContainer* m_truthinfo = nullptr;
};

#endif
