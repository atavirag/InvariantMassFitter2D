// MESelector.cpp
#include "MESelector.h"
#include "MEPair.h"
#include "Math/Vector4D.h"
#include <TFile.h>
#include <TObject.h>
#include <TTree.h>
#include <Rtypes.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <tuple>

constexpr double pdgMass = 1.86484;  // GeV/c^2, PDG D0 mass

// Constructor
MESelector::MESelector(
    const std::tuple<
        std::vector<MECandidates>,
        std::vector<MECandidates>
    >& rawEntries
)
    : _rawEntries(rawEntries)
{
}

void MESelector::applyCuts(double ptMin, double ptMax, double etaMin, double etaMax,
                           double massMin, double massMax, bool removeAmbiguous) {
    std::get<0>(_selectedEntries).clear();
    std::get<1>(_selectedEntries).clear();

    auto select = [&](const std::vector<MECandidates>& raw, std::vector<MECandidates>& sel) {
        for (const auto& c : raw) {
            if (c.pt < ptMin || c.pt > ptMax) continue;
            if (c.eta < etaMin || c.eta > etaMax) continue;
            if (c.m < massMin || c.m > massMax) continue;
            if (removeAmbiguous && TESTBIT(c.type, SelectedD) && TESTBIT(c.type, SelectedDbar)) continue;
            sel.push_back(c);
        }
    };
    select(std::get<0>(_rawEntries), std::get<0>(_selectedEntries));
    select(std::get<1>(_rawEntries), std::get<1>(_selectedEntries));
    std::cout << "Selected " << std::get<0>(_selectedEntries).size() << " / "
              << std::get<1>(_selectedEntries).size() << " candidates (tree 1 / tree 2)\n";
}

// The producer writes each Cand2 once per preceding Cand1 in the collision
void MESelector::deduplicate() {
    auto dedup = [](std::vector<MECandidates>& v) {
        auto tie = [](const MECandidates& c) {
            return std::tie(c.timestamp, c.gIndexCol, c.pt, c.eta, c.phi, c.m, c.type);
        };
        std::sort(v.begin(), v.end(), [&](const auto& a, const auto& b) { return tie(a) < tie(b); });
        v.erase(std::unique(v.begin(), v.end(),
                            [&](const auto& a, const auto& b) { return tie(a) == tie(b); }),
                v.end());
    };
    dedup(std::get<0>(_selectedEntries));
    dedup(std::get<1>(_selectedEntries));
}

void MESelector::buildMixedPairs(int mixDepth) {
    _pairs.clear();
    using EventKey = std::pair<int64_t, int>;
    using EventMap = std::map<EventKey, std::vector<const MECandidates*>>;
    std::map<int, EventMap> c1ByPool, c2ByPool;   // poolBin -> event -> candidates

    for (const auto& c : std::get<0>(_selectedEntries))
        c1ByPool[c.poolBin][{c.timestamp, c.gIndexCol}].push_back(&c);
    for (const auto& c : std::get<1>(_selectedEntries))
        c2ByPool[c.poolBin][{c.timestamp, c.gIndexCol}].push_back(&c);

    auto wrapPhi = [](double dphi) {
        while (dphi < -M_PI / 2) dphi += 2 * M_PI;
        while (dphi >= 3 * M_PI / 2) dphi -= 2 * M_PI;
        return dphi;
    };

    for (const auto& [pool, events1] : c1ByPool) {
        auto it = c2ByPool.find(pool);
        if (it == c2ByPool.end()) continue;

        std::vector<EventKey> keys2;
        std::vector<const std::vector<const MECandidates*>*> cands2;
        for (const auto& [k, v] : it->second) { keys2.push_back(k); cands2.push_back(&v); }

        size_t start = 0;                                   // rotate to spread the mixing
        for (const auto& [k1, v1] : events1) {
            int used = 0;
            for (size_t j = 0; j < keys2.size() && used < mixDepth; ++j) {
                size_t idx = (start + j) % keys2.size();
                if (keys2[idx] == k1) continue;             // never mix within the same event
                for (const auto* a : v1)
                    for (const auto* b : *cands2[idx]) {
                        MEPair p{};
                        p.pt1 = a->pt; p.eta1 = a->eta; p.phi1 = a->phi; p.m1 = a->m; p.type1 = a->type;
                        p.pt2 = b->pt; p.eta2 = b->eta; p.phi2 = b->phi; p.m2 = b->m; p.type2 = b->type;
                        p.poolBin = pool;
                        p.deltaEta = a->eta - b->eta;
                        p.deltaPhi = wrapPhi(a->phi - b->phi);
                        _pairs.push_back(p);
                    }
                ++used;
            }
            ++start;
        }
    }
    std::cout << "Built " << _pairs.size() << " mixed pairs\n";
}

void MESelector::saveOutputTrees(const std::string& outputname) {
    TFile out(outputname.c_str(), "RECREATE");

    if (out.IsZombie()) {
        std::cerr << "Could not open " << outputname << std::endl;
        return;
    }

    TTree tree("pairs", "Mixed-event candidate pairs");

    MEPair p{};

    // Candidate 1
    tree.Branch("pt1",   &p.pt1,   "pt1");
    tree.Branch("eta1",  &p.eta1,  "eta1");
    tree.Branch("phi1",  &p.phi1,  "phi1");
    tree.Branch("m1",    &p.m1,    "m1");
    tree.Branch("y1",    &p.y1,    "y1");
    tree.Branch("type1", &p.type1, "type1/b");

    // Candidate 2
    tree.Branch("pt2",   &p.pt2,   "pt2");
    tree.Branch("eta2",  &p.eta2,  "eta2");
    tree.Branch("phi2",  &p.phi2,  "phi2");
    tree.Branch("m2",    &p.m2,    "m2");
    tree.Branch("y2",    &p.y2,    "y2");
    tree.Branch("type2", &p.type2, "type2/b");

    // Pair-level quantities
    double pairMass = 0.0;
    double pairPt = 0.0;

    tree.Branch("pairMass", &pairMass, "pairMass/D");
    tree.Branch("pairPt",   &pairPt,   "pairPt/D");

    // Angular correlations
    tree.Branch("deltaEta", &p.deltaEta, "deltaEta/D");
    tree.Branch("deltaPhi", &p.deltaPhi, "deltaPhi/D");

    // Fill tree
    for (const auto& pair : _pairs) {
        p = pair;

        // Candidate rapidities, calculated with the nominal D0 mass
        ROOT::Math::PtEtaPhiMVector v1(
            p.pt1, p.eta1, p.phi1, pdgMass
        );

        ROOT::Math::PtEtaPhiMVector v2(
            p.pt2, p.eta2, p.phi2, pdgMass
        );

        p.y1 = v1.Rapidity();
        p.y2 = v2.Rapidity();

        // Pair four-momentum
        const auto vPair = v1 + v2;

        pairMass = vPair.M();
        pairPt   = vPair.Pt();

        tree.Fill();
    }

    out.cd();
    tree.Write("", TObject::kOverwrite);
    out.Close();

    std::cout << "Saved " << tree.GetEntries()
              << " mixed-event pairs to " << outputname << std::endl;
}