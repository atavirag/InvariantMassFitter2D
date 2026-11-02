#pragma once
#include <string>
#include <vector>

class RooWorkspace;

struct MESCorrConfig {
    int    nDeta    = 28;
    double detaMax  = 1.4;
    int    nDphi    = 32;
    double normDeta = 0.2;   // plateau window for the ME normalization
    double normDphi = 0.4;
};

class MESWeights {
public:
    using CorrConfig = MESCorrConfig;
    // Loads the m1 and m2 prefit models of one candidate-pT bin from a fit workspace file.
    // Bins must be added in increasing pT order.
    bool addBin(double ptLo, double ptHi, const std::string& fitFile);

    double weight1(double m, double pt) const;   // candidate 1
    double weight2(double m, double pt) const;   // candidate 2
    double pairWeight(double m1, double pt1, double m2, double pt2) const {
        return weight1(m1, pt1) * weight2(m2, pt2);
    }

    int    binOf(double pt) const { return findBin(pt); }
    double signalFraction1(int bin) const;     // nSig / (nSig + nBkg) of the m1 fit in that pT bin
    double signalFraction2(int bin) const;
    bool   inRange1(int bin, double m) const;  // inside the fit window
    bool   inRange2(int bin, double m) const;

    // Copies the ME pair tree and adds sw1, sw2, swPair. Also fills the (deltaEta, deltaPhi)
    // histograms (raw, weighted, normalized), stores them in the output file and, if
    // plotPrefix is not empty, saves a summary plot.
    bool writeWeightedTree(const std::string& inFile, const std::string& treeName,
                           const std::string& outFile,
                           const CorrConfig& cc = CorrConfig(),
                           const std::string& plotPrefix = "") const;

    // sWeight versus mass for every pT bin (m1 left, m2 right)
    void plotWeightFunctions(const std::string& outputPrefix) const;

    int nBins() const { return (int)_bins.size(); }

    struct WeightFunction {
        bool ok = false;
        double mLo = 0., mHi = 0., step = 0.;
        double nSig = 0., nBkg = 0.;
        std::vector<double> w;                 // weight on a uniform mass grid
        double eval(double m) const;           // linear interpolation, 0 outside the fit range
    };
    struct Bin { double ptLo, ptHi; WeightFunction f1, f2; };

private:

    static bool build(RooWorkspace* ws, const std::string& tag, WeightFunction& out);
    int findBin(double pt) const;

    std::vector<Bin> _bins;
};