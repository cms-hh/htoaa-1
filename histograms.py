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
    output['hLeadingFatJetParticleNetMD_Xbb'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(200, 0., 1, name="hLeadingFatJetParticleNetMD_Xbb", label="hLeadingFatJetParticleNetMD_Xbb", overflow=True,underflow=True)
        .Weight()
    )
    output['hLeadingFatJetDeepTagMD_bbvsLight'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(200, 0., 1, name="hLeadingFatJetDeepTagMD_bbvsLight", label="hLeadingFatJetDeepTagMD_bbvsLight", overflow=True,underflow=True)
        .Weight()
    )
    output['hLeadingFatJetPt'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(200, 0, 1000, name="hLeadingFatJetPt", label="hLeadingFatJetPt", overflow=True,underflow=True)
        .Weight()
    )
    output['hLeadingFatJetMass'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(300, 0, 300, name="hLeadingFatJetMass", label="hLeadingFatJetMass", overflow=True,underflow=True)
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
    output['ak4jetsize_notcleaned'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(10, -0.5, 9.5, name="ak4jetsize_notcleaned", label="ak4jetsize", overflow=True,underflow=True)
        .Weight()
    )
    output['ak4jetsize_notcleaned_mediumbtag'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(10, -0.5, 9.5, name="ak4jetsize_notcleaned_mediumbtag", label="ak4jetsize_mediumbtag", overflow=True,underflow=True)
        .Weight()
    )
    output['dr_lep_ak4jet'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, 0., 5., name="dr_lep_ak4jet", label="dr_lep_ak4jet", overflow=True,underflow=True)
        .Weight()
    )
    output['tau3'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, 0., 1., name="tau3", label="tau3", overflow=True,underflow=True)
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
    output['nConstituents'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 100., name="nConstituents", label="nConstituents", overflow=True,underflow=True)
        .Weight()
    )
    output['subjetpt'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 500., name="subjetpt", label="subjetpt", overflow=True,underflow=True)
        .Weight()
    )
    output['subjet2pt'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 500., name="subjet2pt", label="subjet2pt", overflow=True,underflow=True)
        .Weight()
    )
    output['subjet1_btag'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="subjet1_btag", label="subjet1_btag", overflow=True,underflow=True)
        .Weight()
    )
    output['subjet2_btag'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 1., name="subjet2_btag", label="subjet2_btag", overflow=True,underflow=True)
        .Weight()
    )
    output['subjet1_nBhadrons'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(5, 0., 5., name="subjet1_nBhadrons", label="subjet1_nBhadrons", overflow=True,underflow=True)
        .Weight()
    )
    output['subjet2_nBhadrons'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(5, 0., 5., name="subjet2_nBhadrons", label="subjet2_nBhadrons", overflow=True,underflow=True)
        .Weight()
    )
    output['deltar'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(25, 0., 0.85, name="deltar", label="deltar", overflow=True,underflow=True)
        .Weight()
    )
    output['dq_1'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(50, 0, 5, name="dq_1", label="dq_1", overflow=True,underflow=True)
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
    return output
