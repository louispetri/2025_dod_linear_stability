#ifndef DUNE_SUBTRIANGULATION_DOMAIN_HH
#define DUNE_SUBTRIANGULATION_DOMAIN_HH

#include <vector>
#include <algorithm>
#include <dune/subtriangulation/simpletpmctriangulation/interface.hh>

namespace Dune
{
  namespace SubTriangulation
  {
    namespace DomainDetail
    {
      inline int positionToTag(InterfaceRelativePosition pos)
      {
        switch (pos) {
        case InterfaceRelativePosition::interior:
          return -1;
        case InterfaceRelativePosition::exterior:
          return 1;
        default:
          return 0;
        }
      }

      inline bool compatiblePosition(
          InterfaceRelativePosition domainPos, InterfaceRelativePosition cutCellPos)
      {
        return domainPos == InterfaceRelativePosition::any || domainPos == cutCellPos;
      }
    }

    class Domain
    {
    public:
      using DomainTag = std::vector<int>;
      using Index = std::size_t;

      explicit Domain(Index index)
          : index_(index)
      {
      }

      Domain(Index index, const std::vector<std::string>& positions)
          : index_(index)
      {
        for (const auto& p : positions) {
          addPositionForInterfaces(p);
        }
      }

      bool contains(const std::vector<InterfaceRelativePosition>& interfaceRelativePositions) const
      {
        for (const auto& p : positions_) {
          bool compat = std::equal(p.begin(), p.end(),
              interfaceRelativePositions.begin(), DomainDetail::compatiblePosition);
          if (compat)
            return true;
        }
        return false;
      }

      /*void setPositionForInterface(std::size_t interface, InterfaceRelativePosition pos)
      {
        if (interface >= position_.size())
          position_.resize(interface + 1, InterfaceRelativePosition::any);
        position_[interface] = pos;
      }

      void setPositionForInterfaces(const std::string& v)
      {
        for (std::size_t i = 0; i < v.size(); ++i) {
          if (v[i] == 'i') {
            setPositionForInterface(i, InterfaceRelativePosition::interior);
          } else if (v[i] == 'e') {
            setPositionForInterface(i, InterfaceRelativePosition::exterior);
          } else {
            setPositionForInterface(i, InterfaceRelativePosition::any);
          }
        }
      }*/

      void addPositionForInterfaces(const std::string& v)
      {
        std::vector<InterfaceRelativePosition> position;
        for (char c : v) {
          switch (c) {
          case 'i':
            position.push_back(InterfaceRelativePosition::interior);
            break;
          case 'e':
            position.push_back(InterfaceRelativePosition::exterior);
            break;
          default:
            position.push_back(InterfaceRelativePosition::any);
            break;
          }
        }
        positions_.push_back(position);
      }

      int oldTag() const
      {
        return 0;
      }

      DomainTag tag() const
      {
        return std::vector<int>(positions_[0].size(), 0);
      }

      Index index() const
      {
        return index_;
      }

    private:
      Index index_;
      std::vector<std::vector<InterfaceRelativePosition> > positions_;
    };

    inline bool operator==(const Dune::SubTriangulation::Domain& a, const Dune::SubTriangulation::Domain& b)
    {
      return a.index() == b.index();
    }
  }
}
#endif // DUNE_SubTriangulation_DOMAIN_HH
