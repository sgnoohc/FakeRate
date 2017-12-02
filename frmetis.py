#!/bin/env python

import time
import json

from metis.Sample import DirectorySample
from metis.CondorTask import CondorTask

from metis.StatsParser import StatsParser

import sys
import os

# Configurations
job_tag = "fr_ptratio_study_v1"
exec_path = "scripts/run.sh"
tar_path = "package.tar.gz"
hadoop_path = "metis/fr/{}".format(job_tag)
args = "-c ScanChain.C output.root t -1 dummy" # dummy arguments are there because the executable run.sh was copied from another framework.

# Get into the directory where this lepmetis.py sits. So we can tar up the condor package.
os.system("tar -cf package.tar LeptonTree.cc LeptonTree.h ScanChain.C ScanChain.h pu_weights.h rooutil/*.so rooutil/*.h")
os.chdir("scripts")
os.system("tar -rf ../package.tar *.sh *.C ")
os.chdir("../")
os.system("gzip -f package.tar")

dslocs = [

    ["/QCD_Pt_15to20_MuEnrichedPt5"    , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_20to30_MuEnrichedPt5"    , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_30to50_MuEnrichedPt5"    , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_50to80_MuEnrichedPt5"    , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_80to120_MuEnrichedPt5"   , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v3_2016_fakerate_sample_v2"],
    ["/QCD_Pt_120to170_MuEnrichedPt5"  , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_170to300_MuEnrichedPt5"  , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_300to470_MuEnrichedPt5"  , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_470to600_MuEnrichedPt5"  , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_600to800_MuEnrichedPt5"  , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_800to1000_MuEnrichedPt5" , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_1000toInf_MuEnrichedPt5" , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v3_2016_fakerate_sample_v2"],
    ["/QCD_Pt_20to30_EMEnriched"       , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_30to50_EMEnriched"       , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_50to80_EMEnriched"       , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_80to120_EMEnriched"      , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_120to170_EMEnriched"     , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_170to300_EMEnriched"     , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_300toInf_EMEnriched"     , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_15to20_bcToE"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_15to20_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_20to30_bcToE"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_20to30_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_30to80_bcToE"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_30to80_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_80to170_bcToE"           , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_80to170_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_170to250_bcToE"          , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_170to250_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/QCD_Pt_250toInf_bcToE"          , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/QCD_Pt_250toInf_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_2016_fakerate_sample_v2"],
    ["/Run2016B_DoubleEG"              , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016B_DoubleEG_MINIAOD_03Feb2017_ver2-v2_2016_fakerate_sample_v2"],
    ["/Run2016C_DoubleEG"              , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016C_DoubleEG_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016D_DoubleEG"              , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016D_DoubleEG_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016E_DoubleEG"              , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016E_DoubleEG_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016F_DoubleEG"              , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016F_DoubleEG_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016G_DoubleEG"              , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016G_DoubleEG_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016H_DoubleEG"              , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016H_DoubleEG_MINIAOD_03Feb2017_ver2-v1_2016_fakerate_sample_v2"],
    ["/Run2016H_DoubleEG"              , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016H_DoubleEG_MINIAOD_03Feb2017_ver3-v1_2016_fakerate_sample_v2"],
    ["/Run2016B_DoubleMuon"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016B_DoubleMuon_MINIAOD_03Feb2017_ver2-v2_2016_fakerate_sample_v2"],
    ["/Run2016C_DoubleMuon"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016C_DoubleMuon_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016D_DoubleMuon"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016D_DoubleMuon_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016E_DoubleMuon"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016E_DoubleMuon_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016F_DoubleMuon"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016F_DoubleMuon_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016G_DoubleMuon"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016G_DoubleMuon_MINIAOD_03Feb2017-v1_2016_fakerate_sample_v2"],
    ["/Run2016H_DoubleMuon"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016H_DoubleMuon_MINIAOD_03Feb2017_ver2-v1_2016_fakerate_sample_v2"],
    ["/Run2016H_DoubleMuon"            , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/Run2016H_DoubleMuon_MINIAOD_03Feb2017_ver3-v1_2016_fakerate_sample_v2"],
    ["/DY"                             , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/DY_2016_fakerate_sample_v2"],
    ["/WJets"                          , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/WJets_2016_fakerate_sample_v2"],
    ["/TTbar"                          , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/TTbar_2016_fakerate_sample_v2"],
    ["/WW"                             , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/WW_2016_fakerate_sample_v2"],
    ["/WZ"                             , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/WZ_2016_fakerate_sample_v2"],
    ["/ZZ"                             , "/hadoop/cms/store/user/phchang/metis/lepbaby/2016_fakerate_sample_v2/ZZ_2016_fakerate_sample_v2"],

]

total_summary = {}
while True:
    allcomplete = True
    for ds,loc in dslocs:
        task = CondorTask(
                sample = DirectorySample( dataset=ds, location=loc ),
                open_dataset = False,
                flush = True,
                files_per_output = 20,
                output_name = "merged.root",
                tag = job_tag,
                cmssw_version = "CMSSW_9_2_1", # doesn't do anything
                arguments = args,
                executable = exec_path,
                tarfile = tar_path,
                special_dir = hadoop_path
                )
        task.process()
        allcomplete = allcomplete and task.complete()
        # save some information for the dashboard
        total_summary[ds] = task.get_task_summary()
    # parse the total summary and write out the dashboard
    StatsParser(data=total_summary, webdir="~/public_html/dump/frmetis/").do()
    os.system("chmod -R 755 ~/public_html/dump/frmetis")
    if allcomplete:
        print ""
        print "Job={} finished".format(job_tag)
        print ""
        break
    print "Sleeping 30 seconds ..."
    time.sleep(30)

