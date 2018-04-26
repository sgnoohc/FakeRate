// ROOT
#include "Math/VectorUtil.h"
#include "TVector2.h"

// SNT
#include "rooutil/looper.h"
#include "rooutil/autohist.h"
#include "rooutil/ttreex.h"
#include "rooutil/dorky.h"
#include "pu_weights.h"

#include "TString.h"

#include "LeptonTree.cc"

#include <vector>
#include <tuple>

using namespace std;
using namespace LeptonTreeNameSpace;

struct Lepton
{
    int evt_event;
    int evt_lumiBlock;
    int evt_run;
    int nvtx;
    float instantLumi;
    LV p4;
    float eta;
    int id;
    int nFOs_VVV;
    int njets40;
    int njets40_up;
    int njets40_dn;
    int prescale;
    bool pass_trig;
    float weight;
    float weight_raw;
    std::vector<float> wgt_syst;
    float evt_met;
    float evt_metPhi;
    float evt_met_up;
    float evt_metPhi_up;
    float evt_met_dn;
    float evt_metPhi_dn;
    float evt_mt;
    float evt_mt_up;
    float evt_mt_dn;
    bool passId;
    bool passFO;
    float coneptcorr;
    bool pass_e8i;
    bool pass_m8i;
    bool pass_e17i;
    bool pass_m17i;
    int motherID;
    bool isQCD;
    bool isTTbar;
    bool isEWK;
    bool isData;
    bool isDoubleMuon;
    float reliso;
    float ptratio;
    float ip3d;

    void print()
    {
        std::cout << "evt_event;  " << evt_event  << std::endl;
        std::cout << "p4.pt();    " << p4.pt()    << std::endl;
        std::cout << "p4.eta();   " << p4.eta()   << std::endl;
        std::cout << "id;         " << id         << std::endl;
        std::cout << "nFOs_VVV;   " << nFOs_VVV   << std::endl;
        std::cout << "njets40;    " << njets40    << std::endl;
        std::cout << "prescale;   " << prescale   << std::endl;
        std::cout << "pass_trig;  " << pass_trig  << std::endl;
        std::cout << "weight;     " << weight     << std::endl;
        std::cout << "evt_met;    " << evt_met    << std::endl;
        std::cout << "evt_metPhi; " << evt_metPhi << std::endl;
        std::cout << "evt_mt;     " << evt_mt     << std::endl;
        std::cout << "passId;     " << passId     << std::endl;
        std::cout << "passFO;     " << passFO     << std::endl;
        std::cout << "coneptcorr; " << coneptcorr << std::endl;
    }
};

double calculateMt(const LorentzVector p4, double met, double met_phi);
float getPt(float pt, bool extrPtRel = false);
float getEta(float eta, float ht, bool extrPtRel = false);
void fill(RooUtil::AutoHist& hists, TString prefix, TString hname, TString regionname, bool ismu, bool passtight, float var, float wgt, float rewgt, int nbins, float min, float max);
void fill(RooUtil::AutoHist& hists, TString prefix, TString hname, TString regionname, bool ismu, bool passtight, float var, float wgt, float rewgt, int nbins, double* bounds);
void fill(RooUtil::AutoHist& hists, TString prefix, TString hname, TString regionname, bool ismu, bool passtight, float var, float vary, float wgt, float rewgt, int nbins, double* bounds, int nbinsy, double* boundsy);
void fillEventLevelHistogramsSyst(std::vector<Lepton> leptons, RooUtil::AutoHist& hists, int njets, float evt_met, float evt_metPhi, float evt_mt, int iwgt, TString);

#include "nvtxreweight.h"
