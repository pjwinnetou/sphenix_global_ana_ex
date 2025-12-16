#include <iostream>
#include <fstream>
#include "Ana.h"
using namespace std;

  Ana::Ana(const std::string& name, const char* outname, bool _isMC)
: SubsysReco(name)
  , detector("HCALIN")
{
  outfilename = Form("%s",outname);
  isMC = _isMC;
}

Ana::~Ana()
{
    //delete hm;
    delete towerntuple;
}

int Ana::Init(PHCompositeNode*)
{ 
  triggeranalyzer = new TriggerAnalyzer();
  triggeranalyzer->UseEmulator(useEmulator);
  m_vertex_type_sel.clear();
  m_vertex_type_sel.push_back(GlobalVertex::MBD);
  try {
    InitOutputFile();
    InitTree();
    std::ifstream infile("/sphenix/user/samfred/projects/mbdt0/histmaking/MbdOut.corr");
    if(!infile) {
      std::cerr << "Cannot open file\n";
      return 1;
    }

    int run;
    double val;
    while (infile >> run >> val) {
      mbd_t0_corr[run] = val;
    }

  }  
  catch (const std::exception& e){
    std::cerr << "Ana::Init - Exception during init! For " << e.what() << " return -1" << std::endl;
    return -1;
  }
  return 0;
}

void Ana::InitOutputFile(){
  std::cout << "output filename : " << outfilename.c_str() << std::endl;
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  if(!outfile || outfile->IsZombie()){
    throw std::runtime_error("Failed open file");
  }
}

void Ana::InitTree(){
    towerntuple = new TTree("towerntup", "Ntuple");
    towerntuple->Branch("RunNumber",&m_run_number);
    towerntuple->Branch("ScaledTriggerBit",ScaledTriggerBit,"ScaledTriggerBit[64]/O");
    towerntuple->Branch("LiveTriggerBit",LiveTriggerBit,"LiveTriggerBit[64]/O");
    towerntuple->Branch("vz",&vz);

    towerntuple->Branch("mbd_nhits_south",&m_mbd_hits_south);
    towerntuple->Branch("mbd_nhits_north",&m_mbd_hits_north);
    towerntuple->Branch("mbd_mean_time",&m_mbd_mean_time);
    towerntuple->Branch("m_mbd_time_south",&m_mbd_time_south);
    towerntuple->Branch("m_mbd_time_north",&m_mbd_time_north);
    towerntuple->Branch("m_mbd_sumcharge_south",&m_mbd_sumcharge_south);
    towerntuple->Branch("m_mbd_sumcharge_north",&m_mbd_sumcharge_north);

    towerntuple->Branch("m_sepd_energy",m_sepd_energy,"m_sepd_energy[768]/F");
    towerntuple->Branch("m_sepd_arm",m_sepd_arm,"m_sepd_arm[768]/I");
    towerntuple->Branch("m_sepd_phi",m_sepd_phi,"m_sepd_phi[768]/F");
    towerntuple->Branch("m_sepd_radius",m_sepd_radius,"m_sepd_radius[768]/F");

    towerntuple->Branch("m_zdc_energy_south",m_zdc_energy_south, "m_zdc_energy_south[4]/F");
    towerntuple->Branch("m_zdc_energy_north",m_zdc_energy_north, "m_zdc_energy_north[4]/F");
    towerntuple->Branch("m_zdcinfo_sum_energy_south",&m_zdcinfo_sum_energy_south);
    towerntuple->Branch("m_zdcinfo_sum_energy_north",&m_zdcinfo_sum_energy_north);
    towerntuple->Branch("m_zdcinfo_vtx",&m_zdcinfo_vtx);

}

int Ana::process_event(PHCompositeNode* topNode)
{
    if(!topNode){
      std::cerr << "Ana::Init - topnode PHCompositeNode not valid! return -1" << std::endl;
      return -1;
    }
    if(count % 100 == 0) std::cout << "event : " << count  << std::endl;
    count++;
    process_towers(topNode);
    return Fun4AllReturnCodes::EVENT_OK;
}

int Ana::process_towers(PHCompositeNode* topNode)
{
  if(!topNode){
    std::cerr << "Ana::process_towers - topnode not valid! return -1" << std::endl;
    return -1;
  }
  
  //------------------------Global Info ---------------------
  ProcessGlobalEventInfo(topNode);
  if(IsVtxCut){
    if(fabs(vz)>vzcut) return Fun4AllReturnCodes::ABORTEVENT;
  }

  //-----------------------MBD-------------------------------
  MbdOut * mbdout =  static_cast<MbdOut*>(findNode::getClass<MbdOut>(topNode, "MbdOut"));
  if(!mbdout) return Fun4AllReturnCodes::ABORTRUN;
  m_mbd_mean_time = mbdout->get_t0();
  m_mbd_time_south = mbdout->get_time(0);
  m_mbd_time_north = mbdout->get_time(1);
  m_mbd_sumcharge_south = mbdout->get_q(0);
  m_mbd_sumcharge_north = mbdout->get_q(1);
  m_mbd_hits_south = mbdout->get_npmt(0);
  m_mbd_hits_north = mbdout->get_npmt(1);

  //-----------------------ZDC----------------------------
  TowerInfoContainer * zdc_offlinetowers = findNode::getClass<TowerInfoContainer>( topNode, m_zdc_node_name );
  std::vector<double> zdcenergyvector;

  int zdcsize = zdc_offlinetowers->size();
  for (int channel = 0; channel < zdcsize; channel++)
  {
    auto tower = zdc_offlinetowers->get_tower_at_channel(channel);
    float zdc_e = tower->get_energy();
    if (TowerInfoDefs::isZDC(channel))
    {
      zdcenergyvector.push_back(zdc_e);
    }
  }
  m_zdc_energy_south[0] = zdcenergyvector[0];
  m_zdc_energy_south[1] = zdcenergyvector[2];
  m_zdc_energy_south[2] = zdcenergyvector[4];
  m_zdc_energy_south[3] = zdcenergyvector[6];
  m_zdc_energy_north[0] = zdcenergyvector[8];
  m_zdc_energy_north[1] = zdcenergyvector[10];
  m_zdc_energy_north[2] = zdcenergyvector[12];
  m_zdc_energy_north[3] = zdcenergyvector[14];

  Zdcinfo* zdcinfo = findNode::getClass<Zdcinfo>( topNode, m_zdcinfo_node_name );
  m_zdcinfo_sum_energy_south = zdcinfo->get_zdc_energy(0);
  m_zdcinfo_sum_energy_north = zdcinfo->get_zdc_energy(1);
  m_zdcinfo_vtx = zdcinfo->get_zvertex();


  //-------------------sEPD--------------------------------
  TowerInfoContainer* sepd_offlinetowers = findNode::getClass<TowerInfoContainer>(topNode, m_sepd_node_name);
  EpdGeom *epdgeom = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
  int sepd_size = sepd_offlinetowers->size();
  for (int channel = 0; channel < sepd_size;channel++)
  {
    TowerInfo* offlinetower = sepd_offlinetowers->get_tower_at_channel(channel);
    if(!offlinetower->get_isGood()) continue;
    unsigned int key = TowerInfoDefs::encode_epd(channel);
    float energy = offlinetower->get_energy();
    int arm = TowerInfoDefs::get_epd_arm(key);
	  float phi = epdgeom->get_phi(key);
	  float radius = epdgeom->get_r(key);
    m_sepd_energy[channel] = energy;
    m_sepd_arm[channel] = arm;
    m_sepd_phi[channel] = phi;
    m_sepd_radius[channel] = radius;
  }

  towerntuple->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}
  

int Ana::ProcessGlobalEventInfo(PHCompositeNode* topNode){

  vz=-999;
  vx=-999;
  vy=-999;
    
  bool isglbvtx=true;
/*
  MbdVertexMap *mbdvtxmap = findNode::getClass<MbdVertexMap>(topNode,"MbdVertexMap");
  if(!mbdvtxmap || mbdvtxmap->empty()){ 
    isglbvtx=false;
  }
  if(isglbvtx){
    MbdVertex *bvertex= nullptr;
    if (mbdvtxmap)
    {
      for (MbdVertexMap::ConstIter mbditer= mbdvtxmap->begin(); mbditer != mbdvtxmap->end(); ++mbditer)
      {
        bvertex = mbditer->second;
      }
      if(!bvertex){std::cout << "could not find globalvtxmap iter :: set vtx to (-999,-999,-999)" << std::endl;}
      else if(bvertex){
        m_vertex = bvertex->get_z();
      }
    }
  }
*/

  GlobalVertexMap *globalvtxmap = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
  if(!globalvtxmap){
    std::cout << " no node.." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if(globalvtxmap->empty()){ 
    isglbvtx=false;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if(isglbvtx){
    std::vector<GlobalVertex*> vertices = globalvtxmap->get_gvtxs_with_type(m_vertex_type_sel);
    if(!vertices.empty())
    {
      if(vertices.at(0))
      {
        vx= vertices.at(0)->get_x();
        vy= vertices.at(0)->get_y();
        vz= vertices.at(0)->get_z();
      }
      if(vertices.size() >1)
      {
        std::cout << "More than one vertex of selected type!" << std::endl;
      }
    }
  }

  Gl1Packet *gl1_packet = findNode::getClass<Gl1Packet>(topNode, "14001");
  if (gl1_packet)
  {
    uint64_t gl1_scaledtriggervector = gl1_packet->lValue(0, "ScaledVector");
    uint64_t gl1_livetriggervector = gl1_packet->lValue(0, "TriggerVector");

    for (int i = 0; i < 64; i++)
    {
      ScaledTriggerBit[i] = ((gl1_scaledtriggervector & 0x1U) == 0x1U);
      LiveTriggerBit[i] = ((gl1_livetriggervector & 0x1U) == 0x1U);
      gl1_scaledtriggervector = (gl1_scaledtriggervector >> 1U) & 0xffffffffU;
      gl1_livetriggervector = (gl1_livetriggervector >> 1U) & 0xffffffffU;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int Ana::End(PHCompositeNode* /*topNode*/)
{
  outfile->cd();
  towerntuple->Write();
  outfile->Close();
  std::cout << "end..! " << std::endl;
  return 0;
}
