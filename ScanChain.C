//  .
// ..: P. Chang, philip@physics.ucsd.edu

#include "ScanChain.h"

//int nptbins = 7;
//double ptbins[8] = {10., 15., 20., 25., 30., 35., 40., 120.};

int nptbins = 4;
double ptbins[5] = {10., 25., 30., 40., 120.};

int nptbins_fine = 7;
double ptbins_fine[8] = {10., 15., 20., 25., 30., 35., 40., 120.};

int nptbins_coarse = 6;
double ptbins_coarse[7] = {10., 15., 20., 25., 30., 35., 70.};

int netabins_mu = 2;
double etabins_mu[3] = {0., 1.2, 2.4};
int netabins = 3;
double etabins_el[4] = {0., 0.8, 1.479, 2.5};

bool isData = false;

int evt_event_to_print = 2101836250;

// Rebuilding event level variables

int prev_evt_event = 0;
int prev_evt_run = 0;
int prev_evt_lumi = 0;
std::vector<Lepton> leptons;

//========================================================================================
//
//
// Main Looper
//
//
//========================================================================================
//________________________________________________________________________________________
void ScanChain(TChain* chain, TString outputname, TString baseopts, int nEvents = -1)
{
    //------------------------------------------------------------------------------------
    // Identifying which sample I am running on based on the output name.
    // Generally not a good idea to have the output name dictate the behavior of the code.
    // But for simplicity I am treating the "baseopts" as "options".
    //------------------------------------------------------------------------------------
    isData = baseopts.Contains("Double");   // made it global so i can use in other functions. (Not pretty though.)
    bool isDoubleMuon = baseopts.Contains("DoubleMuon") || baseopts.Contains("doublemuon");
    bool isEWK = (baseopts.Contains("WJets") || baseopts.Contains("DY") || baseopts.Contains("TTbar_")
                  || baseopts.Contains("WW_") || baseopts.Contains("WZ_") || baseopts.Contains("ZZ_"));
    bool isTTbar = baseopts.Contains("TTbar_");
    bool isQCD = baseopts.Contains("QCD") || baseopts.Contains("TTbarFake");
    bool noMCMatch = (isData || isEWK);
    bool isVV = baseopts.Contains("WW_") || baseopts.Contains("WZ_") || baseopts.Contains("ZZ_");
    bool doSyst = true;
    if (isQCD)
    { doSyst = false; }
    //------------------------------------------------------------------------------------
    // Constants in the analysis
    //------------------------------------------------------------------------------------
    float lumi = 35.87;
    // Hardcoded prescale values for the trigger.
    // Something went wrong with prescales for data.
//    float e8i = 4742;
//    float e17i = 611;
//    float m8i = 3492;
//    float m17i = 171;
    // 2016_fakerate_sample_v8
    float e8i = 4673;
    float e17i = 605;
    float m8i = 3454;
    float m17i = 169;
//    58580.9700839
//    7911.88489799
//    46040.4965593
//    3139.29431998
    //------------------------------------------------------------------------------------
    // Histograms in the analysis
    //------------------------------------------------------------------------------------
    RooUtil::AutoHist hists;
    //------------------------------------------------------------------------------------
    //
    // Main iLoop
    //
    //------------------------------------------------------------------------------------
    RooUtil::Looper<LeptonTree> looper(chain, &lepton_tree, nEvents);
    int ievent = 0;
    while (looper.nextEvent())
    {
        ievent++;
        //--------------------------------
        // Prescale calculation via Z-peak
        //--------------------------------
        if (passes_VVV_cutbased_tight() && p4().pt() > 25. && tag_p4().pt() > 30.)
        {
            float wgt = isData ? 1 : scale1fb() * lumi;
            if (abs(id()) == 11)
            {
                if (HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30() > 0)
                {
                    hists.fill(dilep_mass(), "mll_e8i", wgt, 80, 0, 200);
                }
                if (HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30() > 0)
                {
                    hists.fill(dilep_mass(), "mll_e17i", wgt, 80, 0, 200);
                }
            }
            else if (abs(id()) == 13)
            {
                if (HLT_Mu8_TrkIsoVVL() > 0)
                { hists.fill(dilep_mass(), "mll_m8i", wgt, 80, 0, 200); }
                if (HLT_Mu17_TrkIsoVVL() > 0)
                {
                    hists.fill(dilep_mass(), "mll_m17i", wgt, 80, 0, 200);
                }
            }
        }
        //----------------------
        // Compute jet variables
        //----------------------
        int njets40 = 0;
        int njets40_up = 0;
        int njets40_dn = 0;
        for (unsigned int i = 0; i < jets().size(); i++)
        {
            if (ROOT::Math::VectorUtil::DeltaR(jets()[i], p4()) >= 1.)
            {
                if (jets()[i].pt() > 40. && fabs(jets()[i].eta()) < 2.4)
                { njets40++; }
                if ((isEWK || isData) && doSyst)
                {
                    if ((jets()[i].pt() * (1 + jets_unc()[i])) > 40. && fabs(jets()[i].eta()) < 2.4)
                    { njets40_up++; }
                    if ((jets()[i].pt() * (1 - jets_unc()[i])) > 40. && fabs(jets()[i].eta()) < 2.4)
                    { njets40_dn++; }
                }
            }
        }
        //--------------------------
        // Compute trigger variables
        //--------------------------
        // Compute prescales for QCD and Data.
        // QCD will take prescale from the LeptonTree directly.
        // For data, the prescale values were messed up.
        // Therefore it will use hardcoded values.
        int prescale = 0;
        bool pass_trig = 0; // logical AND with whether it matched
        if (abs(id()) == 11)
        {
            // Check the trigger for our isolated single lepton trigger fired.
            // The trigger has positive value for lepton leg matched.
            // The HLT_* branch value will be negative for leg not matched.
            // (See setHLTBranch() function in CORE for more info.)
            if (p4().pt() >= 10 && p4().pt() < 25)
            {
                prescale = e8i;
                //prescale = HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30();
                pass_trig = true;
                if (HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30() <= 0)
                {
                    hists.fill(__COUNTER__, "counter",  1.,  20, 0, 20);
                    pass_trig = false;
                }
            }
            else if (p4().pt() >= 25)
            {
                prescale = e17i;
                //prescale = HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30();
                pass_trig = true;
                if (HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30() <= 0)
                {
                    hists.fill(__COUNTER__, "counter",  1.,  20, 0, 20);
                    pass_trig = false;
                }
            }
        }
        else if (abs(id()) == 13)
        {
            // Check the trigger for our isolated single lepton trigger fired.
            // The trigger has positive value for lepton leg matched.
            // The HLT_* branch value will be negative for leg not matched.
            // (See setHLTBranch() function in CORE for more info.)
            if (p4().pt() >= 10 && p4().pt() < 25)
            {
                prescale = m8i;
                //prescale = HLT_Mu8_TrkIsoVVL();
                pass_trig = true;
                if (HLT_Mu8_TrkIsoVVL() <= 0)
                {
                    hists.fill(__COUNTER__, "counter",  1.,  20, 0, 20);
                    pass_trig = false;
                }
            }
            else if (p4().pt() >= 25)
            {
                prescale = m17i;
                //prescale = HLT_Mu17_TrkIsoVVL();
                pass_trig = true;
                if (HLT_Mu17_TrkIsoVVL() <= 0)
                {
                    hists.fill(__COUNTER__, "counter",  1.,  20, 0, 20);
                    pass_trig = false;
                }
            }
//            if (p4().pt() < 20.)
//            {
//                pass_trig = false;
//            }
        }
        //---------------------
        // Compute event weight
        //---------------------
        float weight = isData ? prescale : scale1fb() * lumi * getTruePUw_Moriond(nvtx());
        //-------------------
        // compute MT and MET
        //-------------------
        float evt_met = evt_corrMET();
        float evt_metPhi = evt_corrMETPhi();
        float evt_mt = calculateMt(p4(), evt_met, evt_metPhi);
        float evt_met_up = -999;
        float evt_metPhi_up = -999;
        float evt_mt_up = -999;
        float evt_met_dn = -999;
        float evt_metPhi_dn = -999;
        float evt_mt_dn = -999;
        if (doSyst)
            //if ( false )
        {
            evt_met_up = evt_corrMET_up();
            evt_metPhi_up = evt_corrMETPhi_up();
            evt_mt_up = calculateMt(p4(), evt_met_up, evt_metPhi_up);
            evt_met_dn = evt_corrMET_dn();
            evt_metPhi_dn = evt_corrMETPhi_dn();
            evt_mt_dn = calculateMt(p4(), evt_met_dn, evt_metPhi_dn);
        }
        //------------------------------------------------------------
        // Lepton ID selection (most important selection in fake rate)
        //------------------------------------------------------------
        //bool passId = abs(id()) == 13 ? passes_VVV_cutbased_tight() : passes_VVV_cutbased_tight() && threeChargeAgree();
        bool passId = abs(id()) == 13 ? passes_VVV_cutbased_fo_noiso() && (RelIso03EA() < 0.4) && (ptratio() > 0.9) && (fabs(ip3d()) < 0.015) : passes_VVV_cutbased_tight() && threeChargeAgree(); // PTRATIOMETHOD
        //bool passId = abs( id() ) == 13 ? passes_VVV_cutbased_tight() : passes_VVV_cutbased_tight();
        bool passFO = abs(id()) == 13 ? passes_VVV_cutbased_veto() : passes_VVV_cutbased_fo_noiso();
        float relIso = RelIso03EAv2();
        //v6
        passFO = abs(id()) == 13 ? passes_VVV_cutbased_fo_noiso() : passes_VVV_cutbased_fo_noiso() && threeChargeAgree() && fabs(ip3d()) < 0.015;
        //passFO = abs( id() ) == 13 ? passes_VVV_cutbased_fo_noiso() : passes_VVV_cutbased_fo_noiso() && fabs(ip3d())<0.015;
        if (abs(id()) == 11)
        { passFO = passFO && (relIso < 0.2); }
        else if (abs(id()) == 13)
        { passFO = passFO && (relIso < 0.4) && (fabs(ip3d()) < 0.015); }
        //{ passFO = passFO && (relIso < 0.4); } // PTRATIOMETHOD
        // If electron, there are some trigger cuts necessary
        // to stay within the trigger online lepton ID

        //------------------------------------------------------------
        // Above portion of the messy code is kept as an archive
        //------------------------------------------------------------
        // New muon IDs (http://cern.ch/go/pm6R)
        if (abs(id()) == 13)
        {
            passId = passes_VVV_cutbased_fo_noiso()
                     && fabs(p4().eta()) < 2.4
                     && fabs(dxyPV()) <= 0.05
                     && fabs(dZ()) <= 0.1
                     //&& ptratio() > 0.9
                     && ptratio() > 0.855
                     && fabs(ip3d()) < 0.015;
            passFO = passes_VVV_cutbased_fo_noiso()
                     && fabs(p4().eta()) < 2.4
                     && fabs(dxyPV()) <= 0.05
                     && fabs(dZ()) <= 0.1
                     && ptratio() > 0.65
                     && fabs(ip3d()) < 0.015;
        }
        // New electron IDs (http://cern.ch/go/668q)
        // 
        // To write down the WPs here (also in the slide): Proposing (including your input):
        //     for SS:
        //     Barrel: MVA > 0.941, Irel0.4 < 0.05, IP3D < 0.010
        //     Endcap: MVA > 0.925, Irel0.4 < 0.07, IP3D < 0.010
        //     for 3l:
        //     Barrel: MVA > 0.920, Irel0.4 < 0.10, IP3D < 0.015
        //     Endcap: MVA > 0.880, Irel0.4 < 0.10, IP3D < 0.015
        // 
        if (abs(id()) == 11)
        {
            bool isEB = fabs(etaSC()) > 1.479 ? false : true;
            bool presel = passes_VVV_cutbased_veto_noiso_noip();
            //bool presel = passes_VVV_cutbased_fo_noiso() && fabs(dxyPV()) <= 0.05 && fabs(dZ()) <= 0.1 && fabs(ip3d()) / ip3derr() < 4;
            passId = presel;
            passFO = presel;

//            // SS tight
//            if (isEB)
//            {
//                if (!(mva_25ns()     > 0.941)) passId = false;
//                if (!(RelIso04EAv2() < 0.05 )) passId = false;
//                if (!(fabs(ip3d())   < 0.01 )) passId = false;
//                if (!(threeChargeAgree()    )) passId = false;
//            }
//            else
//            {
//                if (!(mva_25ns()     > 0.925)) passId = false;
//                if (!(RelIso04EAv2() < 0.07 )) passId = false;
//                if (!(fabs(ip3d())   < 0.01 )) passId = false;
//                if (!(threeChargeAgree()    )) passId = false;
//            }
//            if (isEB)
//            {
//                if (!(mva_25ns()     > 0.941)) passFO = false;
//                if (!(RelIso04EAv2() < 0.4  )) passFO = false;
//                if (!(fabs(ip3d())   < 0.01 )) passFO = false;
//                if (!(threeChargeAgree()    )) passFO = false;
//            }
//            else
//            {
//                if (!(mva_25ns()     > 0.925)) passFO = false;
//                if (!(RelIso04EAv2() < 0.4  )) passFO = false;
//                if (!(fabs(ip3d())   < 0.01 )) passFO = false;
//                if (!(threeChargeAgree()    )) passFO = false;
//            }

            // 3L loose
            if (isEB)
            {
                if (!(mva_25ns()     > 0.920)) passId = false;
                if (!(RelIso04EAv2() < 0.1  )) passId = false;
                if (!(fabs(ip3d())   < 0.015)) passId = false;
            }
            else
            {
                if (!(mva_25ns()     > 0.880)) passId = false;
                if (!(RelIso04EAv2() < 0.1  )) passId = false;
                if (!(fabs(ip3d())   < 0.015)) passId = false;
            }
            if (isEB)
            {
                if (!(mva_25ns()     > 0.920)) passFO = false;
                if (!(RelIso04EAv2() < 0.4  )) passFO = false;
                if (!(fabs(ip3d())   < 0.01 )) passFO = false;
            }
            else
            {
                if (!(mva_25ns()     > 0.880)) passFO = false;
                if (!(RelIso04EAv2() < 0.4  )) passFO = false;
                if (!(fabs(ip3d())   < 0.01 )) passFO = false;
            }

        }

        if (abs(id()) == 11)
        {
            if (passFO || passId)
            //if (passFO)
            {
                bool isEB = fabs(etaSC()) > 1.479 ? false : true;
                float sIeIe = sigmaIEtaIEta_full5x5();
                float hoe = hOverE();
                float deta = fabs(dEtaIn());
                float dphi = fabs(dPhiIn());
                float invep = fabs(1. / ecalEnergy() - 1. / p4().P());
                float cut_sIeIe = isEB ? 0.011 : 0.031;
                float cut_hoe   = 0.08;
                float cut_deta  = 0.01;
                float cut_dphi  = isEB ? 0.04 : 0.08;
                float cut_invep = 0.01;
                bool passHltCuts = (sIeIe < cut_sIeIe && hoe < cut_hoe && deta < cut_deta
                                    && dphi < cut_dphi && invep < cut_invep);
                float ePFIso = ecalPFClusterIso() / p4().pt();
                float hPFIso = hcalPFClusterIso() / p4().pt();
                float trkIso = tkIso() / p4().pt();
                float cut_ePFIso = 0.45;
                float cut_hPFIso = 0.25;
                float cut_trkIso  = 0.2;
                passHltCuts = passHltCuts && ePFIso < cut_ePFIso && hPFIso < cut_hPFIso
                              && trkIso < cut_trkIso;
                passFO = passHltCuts && passFO;
                passId = passHltCuts && passId;
            }
//            if (passId)
//            {
//                bool isEB = fabs(etaSC()) > 1.479 ? false : true;
//                float sIeIe = sigmaIEtaIEta_full5x5();
//                float hoe = hOverE();
//                float deta = fabs(dEtaIn());
//                float dphi = fabs(dPhiIn());
//                float invep = fabs(1. / ecalEnergy() - 1. / p4().P());
//                float cut_sIeIe = isEB ? 0.011 : 0.031;
//                float cut_hoe   = 0.06;
//                float cut_deta  = isEB ? 0.004 : 999;
//                float cut_dphi  = isEB ? 0.02 : 999;
//                float cut_invep = 0.01;
//                bool passHltCuts = (sIeIe < cut_sIeIe && hoe < cut_hoe && deta < cut_deta
//                                    && dphi < cut_dphi && invep < cut_invep);
//                float ePFIso = ecalPFClusterIso() / p4().pt();
//                float hPFIso = hcalPFClusterIso() / p4().pt();
//                float trkIso = tkIso() / p4().pt();
//                float cut_ePFIso = isEB ? 0.16 : 0.12;
//                float cut_hPFIso = isEB ? 0.12 : 0.12;
//                float cut_trkIso = isEB ? 0.08 : 0.08;
//                passHltCuts = passHltCuts && ePFIso < cut_ePFIso && hPFIso < cut_hPFIso
//                              && trkIso < cut_trkIso;
//                passId = passHltCuts && passId;
//            }
        }
        //---------------------------------------------------------
        // Compute cone correction variable for VVV tight isolation
        //---------------------------------------------------------
        float coneptcorr = 0;
        if (abs(id()) == 11)
        {
//            if (abs(etaSC()) <= 1.479)
//            { coneptcorr = std::max(0., relIso - 0.0588); }
//            else
//            { coneptcorr = std::max(0., relIso - 0.0571); }
            //if (abs(etaSC()) <= 1.479)
            //{ coneptcorr = std::max(0., RelIso04EAv2() - 0.05); }
            //else
            //{ coneptcorr = std::max(0., RelIso04EAv2() - 0.07); }

            coneptcorr = std::max(0., RelIso04EAv2() - 0.10);
        }
        if (abs(id()) == 13)
        //{ coneptcorr = std::max(0., relIso - 0.06); }
        //{ coneptcorr = max(double(0.), ((0.9/ptratio()) - 1.)); } // PTRATIOMETHOD
        { coneptcorr = max(double(0.), ((0.855/ptratio()) - 1.)); } // PTRATIOMETHOD
        // If the lepton passed the tight ID no correction is needed.
        if (passId)
        { coneptcorr = 0; }
        //==========================================================================================
        //
        //
        // Done computing variables
        //
        //
        //==========================================================================================
        //--------------------------
        // Instantiate Lepton object
        //--------------------------
        Lepton lepton;
        lepton.evt_event = lepton_tree.evt_event();
        lepton.evt_lumiBlock = lepton_tree.evt_lumiBlock();
        lepton.evt_run = lepton_tree.evt_run();
        lepton.instantLumi = isQCD ? -999 : lepton_tree.instantLumi();
        lepton.nvtx = lepton_tree.nvtx();
        lepton.p4 = p4();
        lepton.id = id();
        lepton.njets40 = njets40;
        lepton.njets40_up = njets40_up;
        lepton.njets40_dn = njets40_dn;
        lepton.prescale = prescale;
        lepton.pass_trig = pass_trig;
        lepton.weight = weight;
        lepton.weight_raw = isData ? 1 : scale1fb() * lumi * getTruePUw_Moriond(nvtx());
        lepton.evt_met = evt_met;
        lepton.evt_metPhi = evt_metPhi;
        lepton.evt_mt = evt_mt;
        lepton.evt_met_up = evt_met_up;
        lepton.evt_metPhi_up = evt_metPhi_up;
        lepton.evt_mt_up = evt_mt_up;
        lepton.evt_met_dn = evt_met_dn;
        lepton.evt_metPhi_dn = evt_metPhi_dn;
        lepton.evt_mt_dn = evt_mt_dn;
        lepton.passId = passId;
        lepton.passFO = passFO;
        lepton.coneptcorr = coneptcorr;
        lepton.pass_e8i = HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30();
        lepton.pass_e17i = HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30();
        lepton.pass_m8i = HLT_Mu8_TrkIsoVVL();
        lepton.pass_m17i = HLT_Mu17_TrkIsoVVL();
        lepton.motherID = motherID();
        lepton.isQCD = isQCD;
        lepton.isTTbar = isTTbar;
        lepton.isEWK = isEWK;
        lepton.isData = isData;
        lepton.isDoubleMuon = isDoubleMuon;
        lepton.reliso = RelIso03EA();
        lepton.ptratio = ptratio();
        lepton.ip3d = ip3d();
        std::vector<float> wgt_syst;
        lepton.wgt_syst = wgt_syst;
        //-------------------------------
        // Debug print relevant variables
        //-------------------------------
        if (lepton_tree.evt_event() == evt_event_to_print)
        {
            std::cout << "Lepton in this event" << std::endl;
            lepton.print();
        }
        //--------------------------------------------------------------------
        // Rebuiliding data at the per event level instead of per lepton level
        //--------------------------------------------------------------------
        if (
                (prev_evt_event != lepton_tree.evt_event()
             || prev_evt_run != lepton_tree.evt_run()
             || prev_evt_lumi != lepton_tree.evt_lumiBlock())
             && ievent != 1
           )
        {
            if (leptons.size() > 0)
            {
                // If it is a new event, it means it's time to take care of previous event data
                // Do something with it.
                fillEventLevelHistogramsSyst(leptons, hists, leptons[0].njets40, leptons[0].evt_met, leptons[0].evt_metPhi, leptons[0].evt_mt, 0, "");
                fillEventLevelHistogramsSyst(leptons, hists, leptons[0].njets40_up, leptons[0].evt_met_up, leptons[0].evt_metPhi_up, leptons[0].evt_mt_up, 13, "syst13_");
                fillEventLevelHistogramsSyst(leptons, hists, leptons[0].njets40_dn, leptons[0].evt_met_dn, leptons[0].evt_metPhi_dn, leptons[0].evt_mt_dn, 14, "syst14_");
                // After doing something with it, clear it out.
                leptons.clear();
            }
        }
        //--------------------------------------------------------------------
        // saving previous event run/lumi/evt
        //--------------------------------------------------------------------
        if (ievent == 1)
        {
            prev_evt_event = lepton_tree.evt_event();
            prev_evt_run = lepton_tree.evt_run();
            prev_evt_lumi = lepton_tree.evt_lumiBlock();
        }
        else if ( ievent >= 2 && (prev_evt_event != lepton_tree.evt_event() || prev_evt_run != lepton_tree.evt_run() || prev_evt_lumi != lepton_tree.evt_lumiBlock()) )
        {
            prev_evt_event = lepton_tree.evt_event();
            prev_evt_run = lepton_tree.evt_run();
            prev_evt_lumi = lepton_tree.evt_lumiBlock();
        }
        // add the lepton to the leptons vector
        if (passFO)
        { leptons.push_back(lepton); }
        //------------------
        // Print lepton info
        //------------------
        if (lepton_tree.evt_event() == evt_event_to_print)
        {
            std::cout << "=====meas====" << std::endl;
            std::cout << "printing event info for evt_event = " << lepton_tree.evt_event() << std::endl;
            std::cout << " passId " <<  passId << std::endl;
            std::cout << " passFO " <<  passFO << std::endl;
            std::cout << " pt " <<  p4().pt() << std::endl;
            std::cout << " nFOs_SS " <<  nFOs_SS() << std::endl;
            std::cout << " evt_met " <<  evt_met << std::endl;
            std::cout << " evt_mt " <<  evt_mt << std::endl;
            std::cout << " prescale " << prescale << std::endl;
            std::cout << " pass_trig " << pass_trig << std::endl;
            std::cout << " njets40 " << njets40 << std::endl;
            std::cout << "=====" << std::endl;
        }
    }
    TH1F* h_Nt_e  = (TH1F*) hists.get("evt_lvl_histo_ptvarbin_meas_el");
    TH1F* h_Nl_e  = (TH1F*) hists.get("evt_lvl_histo_ptvarbin_loose_meas_el");
    TH1F* h_Nt_mu = (TH1F*) hists.get("evt_lvl_histo_ptvarbin_meas_mu");
    TH1F* h_Nl_mu = (TH1F*) hists.get("evt_lvl_histo_ptvarbin_loose_meas_mu");
    if (h_Nt_e) { h_Nt_e  -> Rebin(h_Nt_e ->GetNbinsX()); }
    if (h_Nl_e) { h_Nl_e  -> Rebin(h_Nl_e ->GetNbinsX()); }
    if (h_Nt_mu) { h_Nt_mu -> Rebin(h_Nt_mu->GetNbinsX()); }
    if (h_Nl_mu) { h_Nl_mu -> Rebin(h_Nl_mu->GetNbinsX()); }
    float Nt_e  = h_Nt_e  ? h_Nt_e  -> GetBinContent(1) : 0.;
    float Nl_e  = h_Nl_e  ? h_Nl_e  -> GetBinContent(1) : 0.;
    float Nt_mu = h_Nt_mu ? h_Nt_mu -> GetBinContent(1) : 0.;
    float Nl_mu = h_Nl_mu ? h_Nl_mu -> GetBinContent(1) : 0.;
    float e_e = Nt_e / (Nl_e);
    float e_mu = Nt_mu / (Nl_mu);
    cout << "\nReco (el): " << "Nt = " << Nt_e << ", Nl = " << Nl_e << ", e = " << e_e << endl;
    cout << "\nReco (mu): " << "Nt = " << Nt_mu << ", Nl = " << Nl_mu << ", e = " << e_mu <<
         endl;
    hists.save(outputname);
}

//________________________________________________________________________________________
void fillEventLevelHistogramsSyst_v1(
    std::vector<Lepton> leptons,
    RooUtil::AutoHist& hists,
    int njets,
    float evt_met,
    float evt_metPhi,
    float evt_mt,
    int iwgt,
    TString prefix
)
{
    // Count number of loose and tight leptons
    int nFOs = 0;
    int nIds = 0;
    for (auto& lepton : leptons)
    {
        if (lepton.passFO) { nFOs++; }
        if (lepton.passId) { nIds++; }
    }

    // If there are not exactly one loose lepton skip
    if (nFOs != 1) return;
}

//________________________________________________________________________________________
float computePtRel(LorentzVector lepp4, LorentzVector jetp4, bool subtractLep)
{
    if (jetp4.pt() == 0)
    { return 0.; }
    if (subtractLep)
    { jetp4 -= lepp4; }
    float dot = lepp4.Vect().Dot(jetp4.Vect());
    float ptrel = lepp4.P2() - dot * dot / jetp4.P2();
    ptrel = ptrel > 0 ? sqrt(ptrel) : 0.0;
    return ptrel;
}

//________________________________________________________________________________________
double calculateMt(const LorentzVector p4, double met, double met_phi)
{
    float phi1 = p4.Phi();
    float phi2 = met_phi;
    float Et1  = p4.Et();
    float Et2  = met;
    return sqrt(2 * Et1 * Et2 * (1.0 - cos(phi1 - phi2)));
}

//________________________________________________________________________________________
float getPt(float pt, bool extrPtRel)
{
    if (!extrPtRel && pt >= 70.)
    { return 69.; }
    if (extrPtRel && pt >= 150.)
    { return 149.; }
    if (pt < 10.)
    {
        return 11.;    //use this if lower FR histo bound is 10.
    }
    return pt;
}

//________________________________________________________________________________________
float getEta(float eta, float ht, bool extrPtRel)
{
    if (extrPtRel)
    {
        if (ht >= 800)
        { return 799; }
        return ht;
    }
    if (fabs(eta) >= 2.4)
    { return 2.3; }
    return fabs(eta);
}

//________________________________________________________________________________________
bool passIsolatedFO(int id, float eta, float disc, float pt)
{
    if (abs(id) == 13)
    { return true; }
    float aeta = fabs(eta);
    // if (aeta < 0.8) return mva > -0.155;
    // if ((aeta >= 0.8 && aeta <= 1.479)) return mva > -0.56;
    // if (aeta > 1.479) return mva > -0.76;
    if (aeta < 0.8)
    { return disc > mvacut(-0.86, -0.96, -0.3, pt); }
    if ((aeta >= 0.8 && aeta <= 1.479))
    { return disc > mvacut(-0.85, -0.96, -0.36, pt); }
    if (aeta > 1.479)
    { return disc > mvacut(-0.81, -0.95, -0.63, pt); }
    return false;
}

//________________________________________________________________________________________
// returns A if pt<ptmin, B if pt>ptmax, and linear interpolation between. if pt<10, use C
float mvacut(float A, float B, float C, float pt_)
{
    float ptmin = 15;
    float ptmax = 25;
    return pt_ > 10 ? std::min(A, std::max(B,
                                           A + (B - A) / (ptmax - ptmin) * (pt_ - ptmin))) : C;
}

