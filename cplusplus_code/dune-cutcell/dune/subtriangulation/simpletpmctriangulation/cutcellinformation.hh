#ifndef DUNE_SUBTRIANGULATION_CUTCELLINFORMATION_HH
#define DUNE_SUBTRIANGULATION_CUTCELLINFORMATION_HH

#include <dune/common/float_cmp.hh>

#include <dune/grid/common/scsgmapper.hh>

#include <dune/subtriangulation/common/cuboid.hh>
#include <dune/subtriangulation/simpletpmctriangulation/domainconfiguration.hh>

namespace Dune
{
  namespace SubTriangulation
  {
    template <class GV, class LGV>
    class CutCellInformation
    {
    public:
      using Grid = typename GV::Traits::Grid;
      using BoundingBox = Dune::GridCuboid<Grid>;
      using Entity = typename BoundingBox::Entity;
      using ctype = typename GV::ctype;
      using Coordinate = Dune::FieldVector<ctype, GV::dimension>;

      struct Information {
        std::size_t domainIndex;
        BoundingBox boundingBox;
        ctype volume;
        bool fillsFundamentalCell;
        std::size_t numberOfSnippets;
        Coordinate centerOfMass;
      };

      CutCellInformation(const GV& fundamentalGridView, const LGV& levelSetGridView,
                         const DomainConfiguration<GV, LGV>& domainConfiguration,
                         unsigned int verbose = 0)
          : gridView_(fundamentalGridView)
          , elementMapper_(gridView_)
          , information_(elementMapper_.size())
          , numberOfCutCells_(domainConfiguration.numberOfDomains())
          , verbose_(verbose)
      {

      }

      template<class LST>
      void init(LST& localSubtriangulation, const DomainConfiguration<GV, LGV>& domainConfiguration)
      {
        std::fill(numberOfCutCells_.begin(), numberOfCutCells_.end(), 0);

        std::size_t count = 0;

        for (const Entity& entity : elements(gridView_)) {
          ++count;
          auto index = elementMapper_.index(entity);
          if (verbose_ > 0 && count % 1000 == 0) {
            std::cout << "\r";
            std::cout << "cutCellInformation: working on element " << index << " / "
                      << gridView_.size(0);
          }

          localSubtriangulation.bindOnVolume(entity);

          for (const Domain& domain : domainConfiguration.domains()) {
            if (localSubtriangulation.hasCutCell(domain.index())) {
              auto rawBoundingBox = localSubtriangulation.boundingBox(domain.index());

              information_[index].emplace_back(Information{ domain.index(),
                BoundingBox(entity, rawBoundingBox.lower,
                  rawBoundingBox.higher, CuboidMode::corners),
                  localSubtriangulation.cellVolume(domain.index()),
                  localSubtriangulation.fillsFundamentalCell(domain.index()),
                  localSubtriangulation.numberOfVolumeSnippets(domain.index()),
                  localSubtriangulation.domainToCenterOfMass(domain.index()) });

              ++numberOfCutCells_[domain.index()];
            }
          }
        }
        if (verbose_ > 0) {
          std::cout << "\n";
        }
      }

      // precondition: cutCellsExist(entity,domainIndex) == true
      const Information& information(const Entity& entity, std::size_t domainIndex) const
      {
        return information(elementMapper_.index(entity), domainIndex);
      }

      std::size_t numberOfDomains(const Entity& entity) const
      {
        return numberOfDomains(elementMapper_.index(entity));
      }

      std::size_t numberOfDomains(std::size_t entityIndex) const
      {
        if (entityIndex >= information_.size()) {
          DUNE_THROW(Dune::Exception, "requested number of domains of entity "
                  << entityIndex << " which has no bounding boxes");
        }
        return information_[entityIndex].size();
      }

      // precondition: cutCellsExist(entity,domainIndex) == true
      const Information& information(std::size_t entityIndex, std::size_t domainIndex) const
      {
        if (entityIndex >= information_.size()) {
          DUNE_THROW(Dune::Exception, "requested bounding box of entity "
                  << entityIndex << " which has no bounding boxes");
        }
        for (const Information& inf : information_[entityIndex]) {
          if (inf.domainIndex == domainIndex)
            return inf;
        }
        DUNE_THROW(Dune::Exception, "requested bounding box of entity "
                << entityIndex << " in domain " << domainIndex << " which does not exist");
      }

      bool filledBySingleCutCell(const Entity& entity) const
      {
        auto index = elementMapper_.index(entity);
        if (index >= information_.size()) {
          return false;
        }
        for (const Information& inf : information_[index]) {
          if (inf.fillsFundamentalCell) {
            return true;
          }
        }
        return false;
      }

      // precondition: filledBySingleCutCell(entity) == true
      std::size_t domainOfSingleCutCell(const Entity& entity) const
      {
        auto index = elementMapper_.index(entity);
        for (const Information& inf : information_[index]) {
          if (inf.fillsFundamentalCell) {
            return inf.domainIndex;
          }
        }
        DUNE_THROW(Dune::Exception, "no single cut cell filling the entity could be found");
      }

      bool cutCellsExist(std::size_t index, std::size_t domainIndex) const
      {
        if (index >= information_.size()) {
          DUNE_THROW(Dune::Exception, "index " << index
                                               << " is not a valid element index (information.size "
                                               << information_.size() << " elementMapper.size "
                                               << elementMapper_.size() << " gridView.size(0) "
                                               << gridView_.size(0));
        }
        for (const Information& inf : information_[index]) {
          if (inf.domainIndex == domainIndex) {
            return true;
          }
        }
        return false;
      }

      bool cutCellsExist(const Entity& entity, std::size_t domainIndex) const
      {
        if (elementMapper_.index(entity) >= elementMapper_.size()) {
          DUNE_THROW(Dune::Exception, "index " << elementMapper_.index(entity)
                                               << " is not a valid element index (information.size "
                                               << information_.size() << " elementMapper.size "
                                               << elementMapper_.size() << " gridView.size(0) "
                                               << gridView_.size(0));
        }
        return cutCellsExist(elementMapper_.index(entity), domainIndex);
      }

      std::size_t numberOfCutCells(std::size_t domainIndex) const
      {
        return numberOfCutCells_[domainIndex];
      }

      void printStatistics() const
      {
        ctype minVolume = std::numeric_limits<ctype>::max();
        ctype minBBVolume = std::numeric_limits<ctype>::max();
        ctype maxVolume = -std::numeric_limits<ctype>::max();
        ctype maxBBVolume = -std::numeric_limits<ctype>::max();
        std::size_t numberFilling = 0;
        std::size_t numberTotal = 0;
        std::size_t numberLarger = 0;
        std::size_t numberWithoutSnippets = 0;
        std::map<std::size_t, unsigned int> domainToCutcells;
        std::map<std::size_t, unsigned int> domainToSnippets;
        for (const auto& infs : information_) {
          for (const Information& inf : infs) {
            ++domainToCutcells[inf.domainIndex];
            domainToSnippets[inf.domainIndex] += inf.numberOfSnippets;
            minVolume = std::min(inf.volume, minVolume);
            maxVolume = std::max(inf.volume, maxVolume);
            auto bbVolume = inf.boundingBox.volume();
            minBBVolume = std::min(bbVolume, minBBVolume);
            maxBBVolume = std::max(bbVolume, maxBBVolume);
            if (Dune::FloatCmp::lt(bbVolume, inf.volume)) {
              ++numberLarger;
              std::cout << "cut cell volume " << inf.volume << " greater than bbox volume " << bbVolume << " ffc = " << inf.fillsFundamentalCell << " diff = " << std::abs(inf.volume-bbVolume) << " ( infs " << infs.size() << " )\n";
            }
            if (inf.numberOfSnippets == 0) {
              ++numberWithoutSnippets;
            }
            numberFilling += inf.fillsFundamentalCell;
            ++numberTotal;
            if (inf.volume == 0.0) {
              std::cout << "got exactly 0 volume but " << inf.numberOfSnippets
                        << " snippets and bbVolume of " << inf.boundingBox.volume() << " in domain "
                        << inf.domainIndex << "\n";
            }
          }
        }
        std::cout << "CutCellInformation: " << numberTotal << " cut cells, " << numberFilling
                  << " fill their fundamental cell\n";
        std::cout << numberLarger << " cutcells have a bigger volume than their bounding boxes\n";
        std::cout << numberWithoutSnippets << " cutcells have no contributing snippets\n";
        std::cout << "minimal volume: " << minVolume << " maximal volume: " << maxVolume
                  << " (ratio: " << maxVolume / minVolume << ")\n";
        std::cout << "minimal bounding box volume: " << minBBVolume
                  << " maximal bounding box volume: " << maxBBVolume
                  << " (ratio: " << maxBBVolume / minBBVolume << ")" << std::endl;
        std::cout << "Number of cutcells per domain:\n";
        for (const auto& v : domainToCutcells) {
          std::cout << "Domain " << v.first << ": " << v.second << " ( "
                    << (100. * v.second) / numberTotal
                    << "% ) number of snippets: " << domainToSnippets[v.first] << std::endl;
        }
      }

    private:
      GV gridView_;
      SingleCodimSingleGeomTypeMapper<GV, 0> elementMapper_;
      std::vector<std::vector<Information> > information_;
      std::vector<std::size_t> numberOfCutCells_;
      unsigned int verbose_;
    };
  }
}

#endif // DUNE_SubTriangulation_CUTCELLINFORMATION_HH
