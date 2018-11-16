import ROOT

myfile = ROOT.TFile.Open("MyTrackShowerFile95.root")
mychain = myfile.Get("MyTrackShowerTree")
nHits = mychain.GetBranch("nHits")
print("nHits")
print("nHits.GetEntries() = " + str(nHits.GetEntries()))

for i in range(nHits.GetEntries()):
	print("nHits.GetEntry(" + str(i) + ") = " + str(nHits.GetEntry(i)))

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("xPositions")
xPositions = mychain.GetBranch("xPositions")
print("xPositions.GetEntries() = " + str(xPositions.GetEntries()))

for i in range(xPositions.GetEntries()):
	print("xPositions.GetEntry(" + str(i) + ") = " + str(xPositions.GetEntry(i)))
