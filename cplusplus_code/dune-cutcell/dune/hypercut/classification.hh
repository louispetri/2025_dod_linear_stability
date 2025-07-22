#ifndef DUNE_HYPERCUT_CLASSIFICATION_HH
#define DUNE_HYPERCUT_CLASSIFICATION_HH

#include <cmath>
#include <vector>

namespace Dune::Hypercut
{

enum class CutCellClass { regular, small };

template <class ST>
class CutCellClassifierInterface
{
public:
    using SubTriangulation = ST;

    virtual void classify(SubTriangulation& subTriangulation) = 0;

    bool isSmallCell(std::size_t elementIndex, int domainIndex) const
    {
        return classes_[domainIndex][elementIndex] == CutCellClass::small;
    }

    virtual ~CutCellClassifierInterface()
    {
    }

protected:
    void insertClass(std::size_t elementIndex, int domainIndex, CutCellClass c)
    {
        classes_[domainIndex][elementIndex] = c;
    }

    std::vector<std::vector<CutCellClass>> classes_;
};

template <class ST>
class TriangleVolumeFractionCutCellClassifier : public CutCellClassifierInterface<ST>
{
    using BaseT = CutCellClassifierInterface<ST>;

public:
    using K = typename BaseT::SubTriangulation::ctype;

    TriangleVolumeFractionCutCellClassifier(K threshold)
        : threshold_(threshold)
    {
    }

    void classify(typename BaseT::SubTriangulation& subTriangulation) override
    {
        minSmallCellFraction_ = 1.0;
        maxSmallCellFraction_ = 0.0;

        const auto& indexSet = subTriangulation.gridView().indexSet();
        this->classes_.resize(subTriangulation.domainConfiguration().numberOfDomains());

        for (int domainIndex = 0; domainIndex < subTriangulation.domainConfiguration().numberOfDomains(); ++domainIndex) {
            this->classes_[domainIndex].resize(subTriangulation.gridView().size(0));
        }

        for (const auto& element : elements(subTriangulation.gridView())) {
            subTriangulation.bindOnVolume(element);
            subTriangulation.createCutCells();

            for (auto cutCellIt = subTriangulation.cutCellsBegin(); cutCellIt != subTriangulation.cutCellsEnd(); ++cutCellIt) {
                // compute volume fraction
                const auto& info = subTriangulation.cutCellInformation().information(element, cutCellIt->domainIndex());
                const auto cutCellVolume = info.volume;
                const auto fundamentalElementVolume = info.boundingBox.entity().geometry().volume();
                const auto volumeFraction = cutCellVolume / fundamentalElementVolume;

                // classify based on threshold
                // only consider triangle cells as small
                if (volumeFraction < threshold_ && cutCellIt->isTriangle()) {
                    this->insertClass(indexSet.index(element), cutCellIt->domainIndex(), CutCellClass::small);
                    minSmallCellFraction_ = std::min(volumeFraction, minSmallCellFraction_);
                    maxSmallCellFraction_ = std::max(volumeFraction, maxSmallCellFraction_);
                } else {
                    this->insertClass(indexSet.index(element), cutCellIt->domainIndex(), CutCellClass::regular);
                }
            }
        }
    }

    K minSmallCellFraction() const
    {
        return minSmallCellFraction_;
    }

    K maxSmallCellFraction() const
    {
        return maxSmallCellFraction_;
    }

private:
    K threshold_;
    K minSmallCellFraction_;
    K maxSmallCellFraction_;
};

template <class ST>
class VolumeFractionCutCellClassifier : public CutCellClassifierInterface<ST>
{
    using BaseT = CutCellClassifierInterface<ST>;

public:
    using K = typename BaseT::SubTriangulation::ctype;

    VolumeFractionCutCellClassifier(K threshold)
        : threshold_(threshold)
    {
    }

    void classify(typename BaseT::SubTriangulation& subTriangulation) override
    {
        const auto& indexSet = subTriangulation.gridView().indexSet();
        this->classes_.resize(subTriangulation.domainConfiguration().numberOfDomains());

        for (int domainIndex = 0; domainIndex < subTriangulation.domainConfiguration().numberOfDomains(); ++domainIndex) {
            this->classes_[domainIndex].resize(subTriangulation.gridView().size(0));
        }

        for (const auto& element : elements(subTriangulation.gridView())) {
            subTriangulation.bindOnVolume(element);
            subTriangulation.createCutCells();

            for (auto cutCellIt = subTriangulation.cutCellsBegin(); cutCellIt != subTriangulation.cutCellsEnd(); ++cutCellIt) {
                const auto& info = subTriangulation.cutCellInformation().information(element, cutCellIt->domainIndex());
                const auto cutCellVolume = info.volume;
                const auto fundamentalElementVolume = info.boundingBox.entity().geometry().volume();
                const auto volumeFraction = cutCellVolume / fundamentalElementVolume;

                if (volumeFraction < threshold_) {
                    this->insertClass(indexSet.index(element), cutCellIt->domainIndex(), CutCellClass::small);
                } else {
                    this->insertClass(indexSet.index(element), cutCellIt->domainIndex(), CutCellClass::regular);
                }
            }
        }
    }

private:
    K threshold_;
};

template <class ST>
class TriangleCartesianFaceFractionCutCellClassifier : public CutCellClassifierInterface<ST>
{
    using BaseT = CutCellClassifierInterface<ST>;

public:
    using K = typename BaseT::SubTriangulation::ctype;

    TriangleCartesianFaceFractionCutCellClassifier(K threshold)
        : threshold_(threshold)
    {
    }

    void classify(typename BaseT::SubTriangulation& subTriangulation) override
    {
        const auto& indexSet = subTriangulation.gridView().indexSet();
        this->classes_.resize(subTriangulation.domainConfiguration().numberOfDomains());

        for (int domainIndex = 0; domainIndex < subTriangulation.domainConfiguration().numberOfDomains(); ++domainIndex) {
            this->classes_[domainIndex].resize(subTriangulation.gridView().size(0));
        }

        for (const auto& element : elements(subTriangulation.gridView())) {
            subTriangulation.bindOnVolume(element);
            subTriangulation.createCutCells();

            for (auto cutCellIt = subTriangulation.cutCellsBegin(); cutCellIt != subTriangulation.cutCellsEnd(); ++cutCellIt) {
                // compute volume fraction
                const auto& info =
                    subTriangulation.cutCellInformation().information(element, cutCellIt->domainIndex());
                const K minCartesianFaceLength = std::max(info.boundingBox.height()[0], info.boundingBox.height()[1]);
                const auto fundamentalElementVolume = info.boundingBox.entity().geometry().volume();
                const auto fundamentalElementFaceLength = std::sqrt(fundamentalElementVolume);
                const auto faceFraction = minCartesianFaceLength / fundamentalElementFaceLength;

                // classify based on threshold
                // only consider triangle cells as small
                if (faceFraction < threshold_ && cutCellIt->isTriangle()) {
                    this->insertClass(indexSet.index(element), cutCellIt->domainIndex(), CutCellClass::small);
                } else {
                    this->insertClass(indexSet.index(element), cutCellIt->domainIndex(), CutCellClass::regular);
                }
            }
        }
    }

private:
    K threshold_;
};

}

#endif