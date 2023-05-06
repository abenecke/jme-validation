from bamboo import treefunctions as op

# Object definitions


def muonDef(mu):
    return op.AND(
        mu.pt >= 5.,
        op.abs(mu.eta) <= 2.4,
        op.abs(mu.dxy) <= 0.05,
        op.abs(mu.dz) <= 0.1,
        mu.miniPFRelIso_all <= 0.4,
        mu.sip3d <= 8,
        mu.looseId
    )


def muonConePt(muons):
    return op.map(muons, lambda lep: op.multiSwitch(
        (op.AND(op.abs(lep.pdgId) != 11, op.abs(lep.pdgId) != 13), lep.pt),
        (op.AND(op.abs(lep.pdgId) == 13, lep.mvaTTH > 0.50), lep.pt),
        0.9*lep.pt*(1.+lep.jetRelIso)
    ))


def elDef(el):
    return op.AND(
        el.pt >= 7.,
        op.abs(el.eta) <= 2.5,
        op.abs(el.dxy) <= 0.05,
        op.abs(el.dz) <= 0.1,
        el.miniPFRelIso_all <= 0.4,
        el.sip3d <= 8,
        # el.mvaNoIso_WPL,
        el.lostHits <= 1
    )


def elConePt(electrons):
    return op.map(electrons, lambda lep: op.multiSwitch(
        (op.AND(op.abs(lep.pdgId) != 11, op.abs(lep.pdgId) != 13), lep.pt),
        # (op.AND(op.abs(lep.pdgId) == 11, lep.mvaTTH > 0.30), lep.pt), # run3 MC don't have mvaTTH for electrons
        (op.AND(op.abs(lep.pdgId) == 11), lep.pt),
        0.9*lep.pt*(1.+lep.jetRelIso)
    ))


def cleanElectrons(electrons, muons):
    cleanedElectrons = op.select(electrons, lambda el: op.NOT(
        op.rng_any(
            muons, lambda mu: op.deltaR(el.p4, mu.p4) <= 0.3))
    )
    return cleanedElectrons




def ak4jetDef(jet):
    return op.AND(
        #jet.jetId & 2,  # tight
        jet.pt > 30.,
        op.abs(jet.eta) <= 2.4
    )

def cleanJets(jets, muons, electrons, sort=True):
    jets = op.select(jets, lambda jet: op.AND(
            op.NOT(op.rng_any(electrons, lambda ele: op.deltaR(jet.p4, ele.p4) < 0.4)),
            op.NOT(op.rng_any(muons, lambda mu: op.deltaR(jet.p4, mu.p4) < 0.4))
        ))
    if sort:
        jets = op.sort(jets, lambda j: -j.pt)
    return jets


def ak8jetDef(jet):
    return op.AND(
        #        jet.jetId & 2,  # tight
        jet.subJet1.isValid,
        jet.subJet2.isValid,
        jet.subJet1.pt > 20.,
        jet.subJet2.pt > 20.,
        op.abs(jet.subJet1.eta) <= 2.4,
        op.abs(jet.subJet2.eta) <= 2.4,
        jet.msoftdrop >= 30.,
        jet.msoftdrop <= 210.,
        jet.pt > 200.,
        op.abs(jet.eta) <= 2.4
    )
