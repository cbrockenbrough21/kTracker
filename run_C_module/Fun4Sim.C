R__LOAD_LIBRARY(libfun4all_vect_file_manager)
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libembedding)
R__LOAD_LIBRARY(libevt_filter)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libSQPrimaryGen)
R__LOAD_LIBRARY(libcalibrator)

int Fun4Sim(const int n_evt = 1000) {
  const double FMAGSTR = -1.044;
  const double KMAGSTR = -1.025;

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 5433);  // Geometry based on run number
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->set_DoubleFlag("SIGX_BEAM", 0.3);
  rc->set_DoubleFlag("SIGY_BEAM", 0.3);
  rc->set_DoubleFlag("Z_UPSTREAM", -700.);

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  // Optional calibration (can be disabled if not needed)
  CalibHitElementPos *cal_ele_pos = new CalibHitElementPos();
  cal_ele_pos->CalibTriggerHit(false);
  se->registerSubsystem(cal_ele_pos);

  // Tracking module with EventReducer
  SQReco *reco = new SQReco();
  reco->Verbosity(100);
  reco->set_legacy_rec_container(true); // keep default true
  reco->set_geom_file_name((string)gSystem->Getenv("E1039_RESOURCE") + "/geometry/geom_run005433.root");
  reco->set_enable_KF(true); // Optional Kalman Filter
  reco->setInputTy(SQReco::E1039);
  reco->setFitterTy(SQReco::KFREF);
  reco->set_evt_reducer_opt("h");  // Enable reducer here - h for hodomask
  reco->set_enable_eval(false);
  reco->set_eval_file_name("eval.root");
  reco->set_enable_eval_dst(false);
  for (int ii = 0; ii <= 3; ii++) reco->add_eval_list(ii);
  se->registerSubsystem(reco);

  // INPUT
  Fun4AllVectEventInputManager *in = new Fun4AllVectEventInputManager("VectIn");
  in->Verbosity(1);
  in->set_tree_name("tree");
  in->fileopen("/project/ptgroup/Catherine/kTracker/data/noisy/MC_JPsi_Pythia8_Target_April17_10000_noisy_onlyElectronic.root");  // <-- Replace with your actual file
  se->registerInputManager(in);

  // OUTPUT
  Fun4AllVectEventOutputManager *out = new Fun4AllVectEventOutputManager("VectOut");
  out->SetFileName("vMC_JPsi_Pythia8_Target_April17_10000_onlyElectronic_Cleaned.root");
  out->SetTreeName("tree"); 
  se->registerOutputManager(out);

  // RUN
  se->run(n_evt);
  se->End();
  se->PrintTimer();

  rc->WriteToFile("recoConsts.tsv");
  std::cout << "All done" << std::endl;

  delete se;
  gSystem->Exit(0);
  return 0;
}
