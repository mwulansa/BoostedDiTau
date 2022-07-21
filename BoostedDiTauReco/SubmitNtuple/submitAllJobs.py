import sys,os

if len(sys.argv) == 1:
    print("Please input either 'first' or 'second'")

iteration = sys.argv[1]

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

if sys.argv[1] == "first" :
    if len(sys.argv) == 2: print("Please specify dataset or background only")
    if len(sys.argv) > 2 and not sys.argv[2].startswith("-"):
        os.system("./makeTarBall.csh")
        dataset = sys.argv[2]
        if dataset == "JetHT": os.system("./submitJobs.csh "+dataset+" -dj") 
        if dataset == "SingleMuon": os.system("./submitJobs.csh "+dataset+" -ds")
        SAMPLE = ['DYJetsToLL_M-50_HTBinned','DYJetsToLL_M-4to50_HTBinned','WJetsToLNu_HTBinned','WW','WZ','ZZ']
        for sample in SAMPLE:
            os.system("./submitJobs.csh "+sample+" -b")
    if "-b" in opts :
        os.system("./makeTarBall.csh")
        SAMPLE = ['DYJetsToLL_M-50_HTBinned','DYJetsToLL_M-4to50_HTBinned','WJetsToLNu_HTBinned','WW','WZ','ZZ']
        for sample in SAMPLE:
            os.system("./submitJobs.csh "+sample+" -b")

if sys.argv[1] == "second" :
    SAMPLE = ['QCD_HTBinned','TTJets']
    for sample in SAMPLE:
        os.system("./submitJobs.csh "+sample+" -b")

