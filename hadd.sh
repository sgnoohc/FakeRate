hadd -f outputs/fakerate_qcd_mu.root $1/QCD_Pt_*Mu*/*root
hadd -f outputs/fakerate_qcd_el.root $1/QCD_Pt_*EM*/*root $1/QCD_Pt_*bcToE*/*root
hadd -f outputs/fakerate_data.root $1/*Double*/*root
hadd -f outputs/fakerate_wj.root $1/*WJ*/*root
hadd -f outputs/fakerate_dy.root $1/*DY*/*root
hadd -f outputs/fakerate_ttbar.root $1/*TTbar*/*root
hadd -f outputs/fakerate_vv.root $1/*WW*/*root $1/*WZ*/*root $1/*ZZ*/*root
