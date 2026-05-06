import uproot

tree = uproot.open("output_ntuples/ttll_601230_mc23a_fullsim.root")["reco"]

for b in tree.keys():
    if "weight" in b:
        print(b)