import ROOT

files = ROOT.std.vector("string")()
files.push_back("output_ntuples/ttll_601230_mc23a_fullsim.root")
files.push_back("output_ntuples/ttll_601230_mc23d_fullsim.root")
files.push_back("output_ntuples/ttll_601230_mc23e_fullsim.root")

df = ROOT.RDataFrame("reco", files)
df.Snapshot("reco", "merged_ttll.root")
