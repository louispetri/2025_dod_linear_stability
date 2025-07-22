#ifndef DUNE_SUBTRIANGULATION_DOMAINCONFIGURATION_HH
#define DUNE_SUBTRIANGULATION_DOMAINCONFIGURATION_HH

#include <cstddef>
#include <vector>
#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/subtriangulation/simpletpmctriangulation/interface.hh>
#include <dune/subtriangulation/simpletpmctriangulation/domain.hh>

namespace Dune
{
  namespace SubTriangulation
  {
    struct NoSuchDomainException : public Dune::Exception
    {
    };

    template <class GV, class LGV>
    class DomainConfiguration
    {
    public:
      using const_interface_iterator = typename std::vector<Interface<LGV> >::const_iterator;
      using const_domain_iterator = typename std::vector<Domain>::const_iterator;

      std::size_t size() const
      {
        return domains_.size();
      }

      void addDomain(const Domain& domain)
      {
        domains_.push_back(domain);
      }

      std::size_t numberOfDomains() const
      {
        return domains_.size();
      }

      std::size_t numberOfInterfaces() const
      {
        return interfaces_.size();
      }

      const_domain_iterator
      findDomain(const std::vector<InterfaceRelativePosition>& interfaceRelativePositions) const
      {
        auto predicate = [&interfaceRelativePositions](const Domain& d) {
          return d.contains(interfaceRelativePositions);
        };

        return std::find_if(domains_.begin(), domains_.end(), predicate);
      }

      const_domain_iterator
      findDomain(std::size_t domainIndex) const
      {
        auto predicate = [&domainIndex](const Domain& d) {
          return d.index() == domainIndex;
        };

        return std::find_if(domains_.begin(), domains_.end(), predicate);
      }

      const_domain_iterator
      findInteriorDomain(std::vector<InterfaceRelativePosition> interfaceRelativePosition) const
      {
        std::replace(interfaceRelativePosition.begin(), interfaceRelativePosition.end(),
                     InterfaceRelativePosition::interface, InterfaceRelativePosition::interior);
        return findDomain(interfaceRelativePosition);
      }

      const_domain_iterator
      findExteriorDomain(std::vector<InterfaceRelativePosition> interfaceRelativePosition) const
      {
        std::replace(interfaceRelativePosition.begin(), interfaceRelativePosition.end(),
                     InterfaceRelativePosition::interface, InterfaceRelativePosition::exterior);
        return findDomain(interfaceRelativePosition);
      }

      void addInterface(const Interface<LGV>& interface)
      {
        interfaces_.push_back(interface);
      }

      const_interface_iterator interfacesBegin() const
      {
        return interfaces_.begin();
      }

      const_interface_iterator interfacesEnd() const
      {
        return interfaces_.end();
      }

      IteratorRange<const_interface_iterator> interfaces() const
      {
        return IteratorRange<const_interface_iterator>(interfacesBegin(), interfacesEnd());
      }

      const_domain_iterator domainsBegin() const
      {
        return domains_.begin();
      }

      const_domain_iterator domainsEnd() const
      {
        return domains_.end();
      }

      IteratorRange<const_domain_iterator> domains() const
      {
        return IteratorRange<const_domain_iterator>(domainsBegin(), domainsEnd());
      }

    private:
      std::vector<Interface<LGV> > interfaces_;
      std::vector<Domain> domains_;
    };
  }
}
#endif // DUNE_SubTriangulation_DOMAINCONFIGURATION_HH
