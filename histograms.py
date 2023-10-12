import hist

def histograms(systList, msoftdrop):
    output = {}
    output['hCutFlow'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(10, 0., 10., name="cutflow", label="cutflow")
        .Weight()
    )
    output['hmet'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="met", label="met")
        .Weight()
    )
    output['hmet_vs_lepmva'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="met", label="met")
        .Reg(100, -1., 1., name="lepmva", label="lepmva")
        .Weight()
    )
    output['hmet_vs_lepiso'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="met", label="met")
        .Reg(100, 0., 5., name="lepiso", label="lepiso")
        .Weight()
    )
    output['hsoftdrop_vs_lepiso'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="softdrop", label="softdrop")
        .Reg(100, 0., 5., name="lepiso", label="lepiso")
        .Weight()
    )
    output['hsoftdrop_vs_lepmva'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(100, 0., 300., name="softdrop", label="softdrop")
        .Reg(100, -1., 1., name="lepmva", label="lepmva")
        .Weight()
    )
    output['hLeadingFatJetMSoftDrop'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        #.Variable([90, 95, 100., 105, 110, 115, 120, 125, 130, 140, 150, 170, 200, 800.], name="hLeadingFatJetMSoftDrop", label="hLeadingFatJetMSoftDrop")
        .Variable([90, 110, 150, 200, 300, 800.], name="hLeadingFatJetMSoftDrop", label="hLeadingFatJetMSoftDrop")
        .Weight()
    )
    output['hLeadingFatJetParticleNetMD_Xbb'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(200, 0., 1, name="hLeadingFatJetParticleNetMD_Xbb", label="hLeadingFatJetParticleNetMD_Xbb")
        .Weight()
    )
    output['hLeadingFatJetDeepTagMD_bbvsLight'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(200, 0., 1, name="hLeadingFatJetDeepTagMD_bbvsLight", label="hLeadingFatJetDeepTagMD_bbvsLight")
        .Weight()
    )
    output['hLeadingFatJetPt'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(200, 0, 1000, name="hLeadingFatJetPt", label="hLeadingFatJetPt")
        .Weight()
    )
    output['hLeadingFatJetMass'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(systList, name="systematic")
        .Reg(300, 0, 300, name="hLeadingFatJetMass", label="hLeadingFatJetMass")
        .Weight()
    )
    output['LHE_HT_gen'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(2500, 0, 3000, name="LHE_HT_gen", label="LHE_HT_gen")
        .Weight()
    )
    output['LHE_Njet_gen'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(5, -0.5, 4.5, name="LHE_Njet_gen", label="LHE_Njet_gen")
        .Weight()
    )
    output['genweight'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(201, -100.5, 100.5, name="genWeight", label="genWeight")
        .Weight()
    )
    output['HT'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Variable([0, 150, 200, 225, 250, 275, 300, 325,  350, 375, 400, 425, 450, 475, 500, 550, 600, 1000.], name="HT", label="HT")
        .Weight()
    )
    output['leptonsize'] = (
        hist.Hist.new
        .StrCat(msoftdrop, name="jet")
        .StrCat(['central'], name="systematic")
        .Reg(2, 0.5, 2.5, name="leptonsize", label="leptonsize")
        .Weight()
    )
    return output
