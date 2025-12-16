#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fstream>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4centrality/PHG4CentralityReco.h>

#include <mbd/MbdReco.h>
#include <globalvertex/GlobalVertexReco.h>
#include <epd/EpdReco.h>
#include <zdcinfo/ZdcReco.h>

#include <centrality/CentralityReco.h>
#include <calotrigger/MinimumBiasClassifier.h>

#include <ana/Ana.h>
#include <sstream>
#include <fstream>
#include <string>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libana.so)
R__LOAD_LIBRARY(libepd.so)
R__LOAD_LIBRARY(libzdcinfo.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4vertex.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libFROG.so)



#endif

void get_scaledowns(int runnumber, int scaledowns[])
{

  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","phnxro","");

  if (db)
  {
    printf("Server info: %s\n", db->ServerInfo());
  }
  else
  {
    printf("bad\n");
  }


  TSQLRow *row;
  TSQLResult *res;
  TString cmd = "";
  char sql[1000];


  for (int is = 0; is < 64; is++)
  {
    sprintf(sql, "select scaledown%02d from gl1_scaledown where runnumber = %d;", is, runnumber);
    printf("%s \n" , sql);

    res = db->Query(sql);

    int nrows = res->GetRowCount();

    int nfields = res->GetFieldCount();
    for (int i = 0; i < nrows; i++) {
      row = res->Next();
      for (int j = 0; j < nfields; j++) {
        scaledowns[is] = stoi(row->GetField(j));
      }
      delete row;
    }

    delete res;
  }
  delete db;
}

void Fun4All_macro_globalexample(const char* infile="DST_CALOFITTING_run2auau_ana509_2024p022_v001-00054543-00000.root", bool isMC=false)
{
    double mbdzvtxcut  = 30;
    string vtxstr = Form("VtxCut%.f",mbdzvtxcut);
    string outdir = "test";
    void * dirf = gSystem->OpenDirectory(outdir.c_str());
    if(dirf) gSystem->FreeDirectory(dirf);
    else {gSystem->mkdir(outdir.c_str(), kTRUE);}
    
    std::cout << "input file : " << infile <<std::endl;
    std::string outputname = std::string("output_") + infile;
    const char *outfile = Form("%s/outtree_%s",outdir.c_str(),outputname.c_str());

    std::cout << "final outfile : " << outfile << std::endl;

    Fun4AllServer *se = Fun4AllServer::instance();
    int verbosity = 0;

    se->Verbosity(verbosity);

    pair<int, int> runseg = Fun4AllUtils::GetRunSegment(infile);
    int runnumber = runseg.first;
    int segment = runseg.second;

    //===============
    // conditions DB flags
    //===============

    // global tag
    recoConsts *rc = recoConsts::instance();
    rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
    rc->set_uint64Flag("TIMESTAMP", runnumber);
    
    //Global Reco  
    MbdReco *mbdreco = new MbdReco();
    se->registerSubsystem(mbdreco);
    EpdReco *epdreco = new EpdReco();
    se -> registerSubsystem(epdreco);
    ZdcReco * zdcreco = new ZdcReco();
    zdcreco -> set_zdc1_cut(0.0);
    zdcreco -> set_zdc2_cut(0.0);
    se -> registerSubsystem( zdcreco );

    GlobalVertexReco *gvertex = new GlobalVertexReco();
    se->registerSubsystem(gvertex);
    
    Fun4AllInputManager *in= new Fun4AllDstInputManager("DST_INPUT");
    in->AddFile(infile);
    se->registerInputManager(in);

    int m_scaledowns[64];
    get_scaledowns(runnumber,m_scaledowns);
    for(int i=0;i<64;i++){
      std::cout << "scaledowns " << i << " : " << m_scaledowns[i] << std::endl;
    }
    
    Ana *ca = new Ana("ana",outfile,isMC);
    ca->SetMbdZVtxCut(mbdzvtxcut);
    ca->SetRunNumber(runnumber);
    ca->SetScaledowns(m_scaledowns);
    se->registerSubsystem(ca);

    std::cout << "now run..." << std::endl;
    se->run(2000);
    se->End();
    std::cout << "ok done.. " << std::endl;

    delete se;
}
