import hist
import numpy as np

def clip_value(values, min_value, max_value):
    values = np.minimum(values, max_value-0.001)
    values = np.maximum(values, min_value)
    return values

def histograms(systList, msoftdrop):
    output = {}
    output['hCutFlow'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(10, 0., 10., name="cutflow", label="cutflow",overflow=True,underflow=True)
        .Weight()
    )
    output['hmet'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="met", label="met", overflow=True,underflow=True)
        .Weight()
    )
    output['hmet_vs_lepmva'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="met", label="met",overflow=True,underflow=True)
        .Reg(100, -1., 1., name="lepmva", label="lepmva", overflow=True,underflow=True)
        .Weight()
    )
    output['hmet_vs_lepiso'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="met", label="met", overflow=True,underflow=True)
        .Reg(100, 0., 5., name="lepiso", label="lepiso", overflow=True,underflow=True)
        .Weight()
    )
    output['hsoftdrop_vs_lepiso'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="softdrop", label="softdrop", overflow=True,underflow=True)
        .Reg(100, 0., 5., name="lepiso", label="lepiso", overflow=True,underflow=True)
        .Weight()
    )
    output['hsoftdrop_vs_lepmva'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="softdrop", label="softdrop", overflow=True,underflow=True)
        .Reg(100, -1., 1., name="lepmva", label="lepmva", overflow=True,underflow=True)
        .Weight()
    )
    output['hLeadingFatJetMSoftDrop'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        #.Variable([90, 95, 100., 105, 110, 115, 120, 125, 130, 140, 150, 170, 200, 800.], name="hLeadingFatJetMSoftDrop", label="hLeadingFatJetMSoftDrop")
        .Reg(100,0.,300., name="hLeadingFatJetMSoftDrop", label="hLeadingFatJetMSoftDrop", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNet_massH_Hto4b_v0'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        #.Variable([90, 95, 100., 105, 110, 115, 120, 125, 130, 140, 150, 170, 200, 800.], name="hLeadingFatJetMSoftDrop", label="hLeadingFatJetMSoftDrop") 
        .Reg(100,0.,300., name="FatJet_particleNet_massH_Hto4b_v0", label="FatJet_particleNet_massH_Hto4b_v0", overflow=True,underflow=True)
        .Weight()
    )
    output['hLeadingFatJetParticleNetMD_Xbb'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(200, 0., 1, name="hLeadingFatJetParticleNetMD_Xbb", label="hLeadingFatJetParticleNetMD_Xbb", overflow=True,underflow=True)
        .Weight()
    )
    output['hLeadingFatJetPt'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(200, 0, 1000, name="hLeadingFatJetPt", label="hLeadingFatJetPt", overflow=True,underflow=True)
        .Weight()
    )
    output['LHE_HT_gen'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(2500, 0, 3000, name="LHE_HT_gen", label="LHE_HT_gen", overflow=True,underflow=True)
        .Weight()
    )
    output['LHE_Njet_gen'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(5, -0.5, 4.5, name="LHE_Njet_gen", label="LHE_Njet_gen", overflow=True,underflow=True)
        .Weight()
    )
    output['genweight'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(201, -100.5, 100.5, name="genWeight", label="genWeight", overflow=True,underflow=True)
        .Weight()
    )
    output['HT'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Variable([0, 150, 200, 225, 250, 275, 300, 325,  350, 375, 400, 425, 450, 475, 500, 550, 600, 1000.], name="HT", label="HT", overflow=True,underflow=True)
        .Weight()
    )
    output['leptonsize'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(2, 0.5, 2.5, name="leptonsize", label="leptonsize", overflow=True,underflow=True)
        .Weight()
    )
    output['ak4jetsize'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(10, -0.5, 9.5, name="ak4jetsize", label="ak4jetsize", overflow=True,underflow=True)
        .Weight()
    )
    output['ak4jetsize_mediumbtag'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(10, -0.5, 9.5, name="ak4jetsize_mediumbtag", label="ak4jetsize_mediumbtag", overflow=True,underflow=True)
        .Weight()
    )
    output['dr_lep_ak4jet'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, 0., 5., name="dr_lep_ak4jet", label="dr_lep_ak4jet", overflow=True,underflow=True)
        .Weight()
    )
    output['tau4_2'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, 0., 1., name="tau4_2", label="tau4_2", overflow=True,underflow=True)
        .Weight()
    )
    output['tau4_3'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, 0., 1., name="tau4_3", label="tau4_3", overflow=True,underflow=True)
        .Weight()
    )
    output['dq_1'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, 0, 5, name="dq_1", label="dq_1", overflow=True,underflow=True)
        .Weight()
    )
    output['dphi_lm'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, -3.14, 3.14, name="dphi_lm", label="dphi_lm", overflow=True,underflow=True)
        .Weight()
    )
    output['dphi_lm_jet'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, -3.14, 3.14, name="dphi_lm_jet", label="dphi_lm_jet", overflow=True,underflow=True)
        .Weight()
    )
    output['sign_deta'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, -3.14, 3.14, name="sign_deta", label="sign_deta", overflow=True,underflow=True)
        .Weight()
    )
    output['db_1'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, 0, 5, name="db_1", label="db_1", overflow=True,underflow=True)
        .Weight()
    )
    output['wq'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(5, 0.5, 5.5, name="wq", label="wq", overflow=True,underflow=True)
        .Weight()
    )
    output['ak4_btag'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, 0, 1., name="ak4_btag", label="ak4_btag", overflow=True,underflow=True)
        .Weight()
    )
    output['ak8jetsize'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(5, 0.5, 5.5, name="ak8jetsize", label="ak8jetsize", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_Haa3b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_Haa3b", label="FatJet_particleNetMD_Hto4b_Haa3b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_Haa4b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_Haa4b", label="FatJet_particleNetMD_Hto4b_Haa4b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_Haa34b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_Haa34b", label="FatJet_particleNetMD_Hto4b_Haa34b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_Haa2b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_Haa2b", label="FatJet_particleNetMD_Hto4b_Haa2b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_Haa01b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_Haa01b", label="FatJet_particleNetMD_Hto4b_Haa01b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_QCD4b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_QCD4b", label="FatJet_particleNetMD_Hto4b_QCD4b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_QCD3b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_QCD3b", label="FatJet_particleNetMD_Hto4b_QCD3b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_QCD2b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_QCD2b", label="FatJet_particleNetMD_Hto4b_QCD2b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_QCD1b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_QCD1b", label="FatJet_particleNetMD_Hto4b_QCD1b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_QCD0b'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_QCD0b", label="FatJet_particleNetMD_Hto4b_QCD0b", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNet_massA_Hto4b_v1'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(200, 0., 80, name="FatJet_particleNet_massA_Hto4b_v1", label="FatJet_particleNet_massA_Hto4b_v1", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNet_massA_Hto4b_v2'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(200, 0., 80., name="FatJet_particleNet_massA_Hto4b_v2", label="FatJet_particleNet_massA_Hto4b_v2", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNet_massA_Hto4b_v3'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(200, 0., 80., name="FatJet_particleNet_massA_Hto4b_v3", label="FatJet_particleNet_massA_Hto4b_v3", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNet_massA_Hto4b_v0'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(200, 0., 80., name="FatJet_particleNet_massA_Hto4b_v0", label="FatJet_particleNet_massA_Hto4b_v0", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNet_massH_Hto4b_v3'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 50., 200., name="FatJet_particleNet_massH_Hto4b_v3", label="FatJet_particleNet_massH_Hto4b_v3", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNet_massH_Hto4b_v2'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 50., 200., name="FatJet_particleNet_massH_Hto4b_v2", label="FatJet_particleNet_massH_Hto4b_v2", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNet_massH_Hto4b_v1'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 50., 200., name="FatJet_particleNet_massH_Hto4b_v1", label="FatJet_particleNet_massH_Hto4b_v1", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNet_massH_Hto4b_v0'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 50., 200., name="FatJet_particleNet_massH_Hto4b_v0", label="FatJet_particleNet_massH_Hto4b_v0", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto4b_vs_QCDsum'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto4b_vs_QCDsum", label="FatJet_particleNetMD_Hto4b_vs_QCDsum", overflow=True,underflow=True)
        .Weight()
    )
    output['FatJet_particleNetMD_Hto34b_vs_QCDsum'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="FatJet_particleNetMD_Hto34b_vs_QCDsum", label="FatJet_particleNetMD_Hto34b_vs_QCDsum", overflow=True,underflow=True)
        .Weight()
    )
    output['apt'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 500., name="pt", label="a_pt", overflow=True,underflow=True)
        .Weight()
    )
    return output
