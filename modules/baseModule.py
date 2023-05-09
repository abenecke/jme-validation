from bamboo.analysismodules import NanoAODModule, HistogramsModule
from bamboo.treedecorators import NanoAODDescription, CalcCollectionsGroups
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo import treefunctions as op
from bamboo.analysisutils import forceDefine

from itertools import chain


class NanoBaseJME(NanoAODModule, HistogramsModule):
    def __init__(self, args):
        super(NanoBaseJME, self).__init__(args)

    def addArgs(self, parser):
        super(NanoBaseJME, self).addArgs(parser)
        parser.add_argument("--era",
                            action='store',
                            type=str,
                            default=None,
                            help='This has no use right now!')

    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        def isMC():
            if sampleCfg['type'] == 'data':
                return False
            elif sampleCfg['type'] == 'mc':
                return True
            else:
                raise RuntimeError(
                    f"The type '{sampleCfg['type']}' of {sample} dataset not understood.")

        era = sampleCfg['era']  # reserved for future use
        self.is_MC = isMC()
        self.triggersPerPrimaryDataset = {}

        def addHLTPath(PD, HLT):
            if PD not in self.triggersPerPrimaryDataset.keys():
                self.triggersPerPrimaryDataset[PD] = []
            try:
                self.triggersPerPrimaryDataset[PD].append(
                    getattr(tree.HLT, HLT))
            except AttributeError:
                print("Couldn't find branch tree.HLT.%s, will omit it!" % HLT)

        def getNanoAODDescription():
            groups = ["HLT_", "MET_","PV_","Pileup_","Rho_"]
            collections = ["nElectron", "nJet", "nMuon", "nFatJet", "nSubJet","nGenJet"]
            varReaders = [CalcCollectionsGroups(Jet=("pt", "mass"))]
            return NanoAODDescription(groups=groups, collections=collections, systVariations=varReaders)

        tree, noSel, backend, lumiArgs = super(NanoBaseJME, self).prepareTree(tree=tree,
                                                                                 sample=sample,
                                                                                 sampleCfg=sampleCfg,
                                                                                 description=getNanoAODDescription(),
                                                                                 backend=backend)
        ### Triggers ###
        # Muon
        addHLTPath('Muon', 'IsoMu24')
        addHLTPath('Muon', 'IsoMu27')
        addHLTPath('Muon', 'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8')
        # EGamma
        addHLTPath('EGamma', 'Ele32_WPTight_Gsf')
        addHLTPath('EGamma', 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL')
        # MuonEG
        addHLTPath('MuonEG', 'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ')
        # SingleMuon
        addHLTPath('SingleMuon', 'IsoMu24')
        addHLTPath('SingleMuon', 'IsoMu27')
        # DoubleMuon
        addHLTPath('DoubleMuon', 'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8')

        # JetHT
        addHLTPath('JetHT', 'DiPFJetAve40')
        addHLTPath('JetHT', 'DiPFJetAve60')
        addHLTPath('JetHT', 'DiPFJetAve80')
        addHLTPath('JetHT', 'DiPFJetAve140')
        addHLTPath('JetHT', 'DiPFJetAve200')
        addHLTPath('JetHT', 'DiPFJetAve260')
        addHLTPath('JetHT', 'DiPFJetAve320')
        addHLTPath('JetHT', 'DiPFJetAve400')
        addHLTPath('JetHT', 'DiPFJetAve500')
        addHLTPath('JetHT', 'DiPFJetAve60_HFJEC')
        addHLTPath('JetHT', 'DiPFJetAve80_HFJEC')
        addHLTPath('JetHT', 'DiPFJetAve100_HFJEC')
        addHLTPath('JetHT', 'DiPFJetAve160_HFJEC')
        addHLTPath('JetHT', 'DiPFJetAve220_HFJEC')
        addHLTPath('JetHT', 'DiPFJetAve300_HFJEC')

        # Gen Weight and Triggers
        if self.is_MC:
            noSel = noSel.refine('genWeight', weight=tree.genWeight, cut=(
                op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values()))))
        else:
            noSel = noSel.refine('trigger', cut=[makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset)])

        #### reapply JECs ###
        from bamboo.analysisutils import configureJets, configureType1MET
        if era == "2022":
            print("isMC ", self.is_MC)
            if self.is_MC: 
                configureJets(tree._Jet, "AK4PFPuppi",
                              jec="Winter22Run3_V2_MC",
                              mayWriteCache=True,
                              isMC=self.is_MC, backend = backend)
                # configureType1MET(tree._MET,
                #     jec="Summer16_07Aug2017_V20_MC",
                #     smear="Summer16_25nsV1_MC",
                #     jesUncertaintySources=["Total"],
                #     mayWriteCache=isNotWorker,
                #     isMC=self.isMC(sample), backend=be)
            else:
                configureJets(tree._Jet, "AK4PFPuppi",
                              jec="Winter22Run3_RunD_V2_DATA",
                              mayWriteCache=True,
                              isMC=self.is_MC, backend = backend)

        for calcProd in tree._Jet.calcProds:
            forceDefine(calcProd,noSel)

        return tree, noSel, backend, lumiArgs
