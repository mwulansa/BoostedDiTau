import os, sys

# assign directory
Sample = sys.argv[1]
#Sample = "TCP_m_30_UL17/Ntuple_TCP_m30_htj_100to400_Summer20UL17_v21"

#Sample = "TCP_m10_ht_400toInf_slc7_amd64_gcc10_MINIAOD"

#directory = '/eos/uscms/store/user/mwulansa/Events/'+Sample
directory = '/eos/uscms/store/user/mwulansa/TCPNtuple/'+Sample+'/'

# iterate over files in
# that directory
for root, dirs, files in os.walk(directory):
    for filename in files:
#        print("root: " + root.replace("/eos/uscms/","root://cmseos.fnal.gov//"))
        print(root.replace("/eos/uscms/","root://cmseos.fnal.gov//")+"/"+filename)
