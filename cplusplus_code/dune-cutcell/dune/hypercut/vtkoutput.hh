#ifndef DUNE_HYPERCUT_VTKOUTPUT_HH
#define DUNE_HYPERCUT_VTKOUTPUT_HH

#include <dune/subtriangulation/io/refinedvtkwriter.hh>
#include <dune/subtriangulation/io/vtkfunction.hh>

namespace Dune::Hypercut {

    template <typename ST>
    class CutCellVolumeUnfittedVTKGridFunction
        : public Dune::SubTriangulation::UnfittedVTKFunction<ST>
    {
    protected:
      using GridView = typename ST::GridView;
      using Entity = typename GridView::template Codim<0>::Entity;
      using Base = Dune::SubTriangulation::UnfittedVTKFunction<ST>;

    public:
      CutCellVolumeUnfittedVTKGridFunction(const ST& _st)
          : st_(_st)
      {
      }

      virtual int ncomps() const
      {
        return 1;
      }

      virtual bool evaluateOn(const typename Base::EntityPartPointer& part) const
      {
        domainIndex_ = part->domainIndex();
        return true;
      }

      virtual bool evaluateOn(const typename Base::IntersectionPartPointer& part) const
      {
        domainIndex_ = part->insideDomainIndex();
        return true;
      }

      virtual double
      evaluate(int comp, const Entity& e,
               const Dune::FieldVector<typename GridView::ctype, GridView::dimension>& xi) const
      {
        if (!st_.cutCellInformation().cutCellsExist(e, domainIndex_))
        {
          return 0.0;
        }

        const auto& info = st_.cutCellInformation().information(e, domainIndex_);
        return info.volume;
      }

      virtual std::string name() const
      {
        return "volume";
      }

    protected:
      const ST& st_;
      mutable int domainIndex_;
    };

    template <class ST, class Classifier>
    class CutCellClassUnfittedVTKGridFunction
        : public Dune::SubTriangulation::UnfittedVTKFunction<ST>
    {
    protected:
      using GridView = typename ST::GridView;
      using Entity = typename GridView::template Codim<0>::Entity;
      using Base = Dune::SubTriangulation::UnfittedVTKFunction<ST>;

    public:
      CutCellClassUnfittedVTKGridFunction(const ST& _st, typename ST::LocalSubTriangulation& localSubTriangulation, Classifier& classifier)
          : st_(_st), localSubTriangulation_(localSubTriangulation), classifier_(classifier)
      {
        classifier_.classify(localSubTriangulation_);
      }

      virtual int ncomps() const
      {
        return 1;
      }

      virtual bool evaluateOn(const typename Base::EntityPartPointer& part) const
      {
        domainIndex_ = part->domainIndex();
        return true;
      }

      virtual bool evaluateOn(const typename Base::IntersectionPartPointer& part) const
      {
        domainIndex_ = part->insideDomainIndex();
        return true;
      }

      virtual double
      evaluate(int comp, const Entity& e,
               const Dune::FieldVector<typename GridView::ctype, GridView::dimension>& xi) const
      {
        if (classifier_.isSmallCell(st_.gridView().indexSet().index(e), domainIndex_))
        {
          return 1.0;
        }

        return 0.0;
      }

      virtual std::string name() const
      {
        return "small cells";
      }

    protected:
      const ST& st_;
      typename ST::LocalSubTriangulation& localSubTriangulation_;
      Classifier& classifier_;
      mutable int domainIndex_;
    };

    template<class ST, class U, class Basis, class LocalVector, int offset, int size>
    class UnfittedVectorVTKGridFunction :
      public Dune::SubTriangulation::UnfittedVTKFunction<ST>
    {
      using Base = Dune::SubTriangulation::UnfittedVTKFunction<ST>;
      using SubTriangulation = typename Base::ST;
      using Entity = typename SubTriangulation::Entity;
      static const int dim = SubTriangulation::dim;
      using DF = typename SubTriangulation::ctype;

    public:
      UnfittedVectorVTKGridFunction(std::string label, const SubTriangulation& subTriangulation, const Basis& basis, const U& u, int domainIndex) :
        label_(label),
        basis_(basis),
        insideView_(basis_.localView(domainIndex)),
        ep_(nullptr),
        ip_(nullptr),
        u_(u),
        domainIndex_(domainIndex),
        outsideDomainIndex_(-1),
        internalIntersection_(false)
      {
        for (int i = 0; i < subTriangulation.domainConfiguration().numberOfDomains(); ++i) {
          outsideViews_.push_back(basis.localView(i));
        }
      }

      virtual bool evaluateOn( const typename Base::EntityPartPointer& part) const
      {
        if (part->domainIndex() != domainIndex_) {
          ep_ = nullptr;
          return false;
        }

        ep_ = part;
        internalIntersection_ = false;
        insideView_.bind(part->entity());
        readLocalVector(u_, insideView_, insideValue_);

        return true;
      }

      virtual bool evaluateOn( const typename Base::IntersectionPartPointer& part) const
      {
        if (part->insideDomainIndex() != domainIndex_ && part->outsideDomainIndex() != domainIndex_) {
          ip_ = nullptr;
          return false;
        }

        ip_ = part;
        internalIntersection_ = false;
        insideView_.bind(part->inside());
        readLocalVector(u_, insideView_, insideValue_);

        internalIntersection_ = false;
        outsideDomainIndex_ = -1;

        if (part->neighbor()) {
          internalIntersection_ = true;
          outsideDomainIndex_ = part->outsideDomainIndex();
          outsideViews_[outsideDomainIndex_].bind(part->outside());
          readLocalVector(u_, outsideViews_[outsideDomainIndex_], outsideValue_);
        }

        return true;
      }

      virtual double evaluate(int comp, const Entity& e, const Dune::FieldVector<DF, dim>& xi) const
      {
        auto boundingBox = ep_ == nullptr ? ip_->boundingBoxInside() : ep_->boundingBox();

        auto insideBasis = insideView_.tree().child(comp + offset).finiteElement().localBasis();
        std::vector<typename std::decay_t<decltype(insideBasis)>::Traits::RangeType> insidePhi(insideBasis.size());
        insideBasis.evaluateFunction(boundingBox.local(e.geometry().global(xi)), insidePhi);

        double y1 = 0.0;

        for (int i = 0; i < insideView_.tree().child(comp + offset).size(); ++i) {
          y1 += insideValue_[insideView_.tree().child(comp + offset).localIndex(i)] * insidePhi[i];
        }

        if (internalIntersection_) {
          auto outsideBasis = outsideViews_[outsideDomainIndex_].tree().child(comp + offset).finiteElement().localBasis();
          std::vector<typename std::decay_t<decltype(outsideBasis)>::Traits::RangeType> outsidePhi(outsideBasis.size());
          outsideBasis.evaluateFunction(ip_->boundingBoxOutside().local(e.geometry().global(xi)), outsidePhi);

          double y2 = 0.0;

          for (int i = 0; i < outsideViews_[outsideDomainIndex_].tree().child(comp + offset).size(); ++i) {
            y2 += outsideValue_[outsideViews_[outsideDomainIndex_].tree().child(comp + offset).localIndex(i)] * outsidePhi[i];
          }

          y1 *= 0.5;
          y1 += 0.5 * y2;
        }

        return y1;
      }

      virtual int ncomps () const
      {
        return size;
      }

      virtual std::string name () const
      {
        return label_;
      }

    private:
      std::string label_;
      const Basis& basis_;
      mutable typename Basis::LocalView insideView_;
      mutable std::vector<typename Basis::LocalView> outsideViews_;
      mutable LocalVector insideValue_;
      mutable LocalVector outsideValue_;
      mutable typename Base::EntityPartPointer ep_;
      mutable typename Base::IntersectionPartPointer ip_;
      const U& u_;
      int domainIndex_;
      mutable int outsideDomainIndex_;
      mutable bool internalIntersection_;
    };

    template<typename ST>
    class UnfittedMultiDomainVTKGridFunction :
      public Dune::SubTriangulation::UnfittedVTKFunction<ST>
    {
      typedef typename ST::GridView GridView;
      typedef typename GridView::ctype DF;
      enum {dim=GridView::dimension};
      typedef typename GridView::template Codim<0>::Entity Entity;

      typedef Dune::SubTriangulation::UnfittedVTKFunction<ST> Base;

    public:

      UnfittedMultiDomainVTKGridFunction
        (const std::string _label) : label_(_label), components_(0), warningPrinted_(false)
      {}

      UnfittedMultiDomainVTKGridFunction* add
        (const Base * func)
      {
        if(components_ != func->ncomps() && components_ > 0)
        {DUNE_THROW(Dune::Exception,"Components of all domain functions must be equal!");}

        components_ = func->ncomps();
        functions_.push_back(BasePointer(func));
        return this;
      }

      virtual int ncomps () const
      {
        return components_;
      }

      virtual bool evaluateOn( const typename Base::EntityPartPointer & part) const
      {
        currentFunction_=-1;
        currentFunctionOut_=-1;
        for(size_t i=0; i<functions_.size(); ++i) {
          if(! functions_[i]->evaluateOn(part) ) continue;
          if(currentFunction_ >= 0)
          {DUNE_THROW(Dune::Exception,"Domains of given function overlap!");}

          currentFunction_ = i;
        }
        return currentFunction_ >= 0;
      }


      virtual bool evaluateOn( const typename Base::IntersectionPartPointer & part) const
      {
        currentFunction_=-1;
        currentFunctionOut_=-1;
        for(size_t i=0; i<functions_.size(); ++i) {
          if(! functions_[i]->evaluateOn(part) ) continue;
          if(currentFunction_ >= 0 && currentFunctionOut_ >= 0)
          {DUNE_THROW(Dune::Exception,"Domains of given function overlap!");}

          currentFunctionOut_ = currentFunction_;
          currentFunction_ = i;

        }
        return currentFunction_ >= 0;
      }

      virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DF,dim>& xi) const
      {
        assert(currentFunction_>=0);
        double y = functions_[currentFunction_]->evaluate(comp,e,xi);
        if(currentFunctionOut_ >= 0)
        {
          if (!warningPrinted_)
          {
            std::cout << "Warning: Evaluating UnfittedMultiDomainVTKGridFunction on interface"
                      << " between two domains - averaging function values..." << std::endl;
            warningPrinted_ = true;
          }
          y *= 0.5;
          y += 0.5*functions_[currentFunctionOut_]->evaluate(comp,e,xi);
        }
        return y;
      }

      virtual std::string name () const
      {
        return label_;
      }

    private:
      const std::string label_;
      typedef std::shared_ptr<const Base> BasePointer;
      std::vector<BasePointer> functions_;
      mutable int currentFunction_;
      mutable int currentFunctionOut_;
      int components_;

      mutable bool warningPrinted_;
    };

      template <typename ST, int m>
      class ExactSolutionUnfittedVTKGridFunction
          : public Dune::SubTriangulation::UnfittedVTKFunction<ST>
      {

        using Base = Dune::SubTriangulation::UnfittedVTKFunction<ST>;
        using SubTriangulation = typename Base::ST;
        using Entity = typename SubTriangulation::Entity;
        static const int dim = SubTriangulation::dim;
        using DF = typename SubTriangulation::ctype;
        using GV = typename ST::GridView;
        enum { n = dim };

      public:
        explicit ExactSolutionUnfittedVTKGridFunction(std::function<Dune::FieldVector<double, m>(Dune::FieldVector<double, n>, double)> sol, DF t)
            : sol_(sol), t_(t)
        {
        }

        virtual int ncomps() const
        {
          return m;
        }

        virtual bool evaluateOn(const typename Base::EntityPartPointer& part) const
        {
          return true;
        }

        virtual bool evaluateOn(const typename Base::IntersectionPartPointer& part) const
        {
          return true;
        }

        virtual double evaluate(int comp, const Entity& e, const Dune::FieldVector<DF, n>& xi) const
        {
          return sol_((e.geometry().global(xi)), t_)[comp];
        }

        virtual std::string name() const
        {
          return "ExactSolution";
        }

      private:
        std::function<Dune::FieldVector<double, m>(Dune::FieldVector<double, n>, double)> sol_;
        DF t_;
      };

template<class ST>
void writeGrid(const ST& subTriangulation, const std::string& filename)
{
    Dune::RefinedVtkWriter<ST, double> vtk(subTriangulation.gridView(), subTriangulation);

    vtk.addVertexData(std::make_shared<Dune::Hypercut::CutCellVolumeUnfittedVTKGridFunction<ST>>(subTriangulation));

    vtk.write(filename, Dune::VTK::OutputType::appendedraw, 2,Dune::UDGVTKWriteMode::writeVolume);
}

template<class ST, class CutCellClassifier>
void writeGrid(const ST& subTriangulation, typename ST::LocalSubTriangulation& localSubTriangulation, CutCellClassifier& classifier, const std::string& filename)
{
    Dune::RefinedVtkWriter<ST, double> vtk(subTriangulation.gridView(), subTriangulation);

    vtk.addVertexData(std::make_shared<Dune::Hypercut::CutCellVolumeUnfittedVTKGridFunction<ST>>(subTriangulation));
    vtk.addVertexData(std::make_shared<Dune::Hypercut::CutCellClassUnfittedVTKGridFunction<ST, CutCellClassifier>>(subTriangulation, localSubTriangulation, classifier));

    vtk.write(filename, Dune::VTK::OutputType::appendedraw, 2,Dune::UDGVTKWriteMode::writeVolume);
}

template<class ST, class Vector, class Basis, int m, int finiteElementSize>
void writeMultiDomainFunction(const ST& subTriangulation, const Vector& u, const Basis& basis, const std::string& filename,
                              std::integral_constant<int, m>, std::integral_constant<int, finiteElementSize>)
{
    using LocalVector = Dune::FieldVector<double, m * finiteElementSize>;

    Dune::RefinedVtkWriter<ST, double> vtk(subTriangulation.gridView(), subTriangulation);

    using UMVGF = UnfittedMultiDomainVTKGridFunction<ST>;
    auto rho = std::make_shared<UMVGF>("rho");
    rho->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 0, 1>(
      "rho0", subTriangulation, basis, u, 0));
    rho->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 0, 1>(
      "rho1", subTriangulation, basis, u, 1));
    vtk.addVertexData(rho);

    if constexpr (m > 1) {
      auto momentum = std::make_shared<UMVGF>("m");
      momentum->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 1, ST::dim>(
        "m0", subTriangulation, basis, u, 0));
      momentum->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 1, ST::dim>(
        "m1", subTriangulation, basis, u, 1));
      vtk.addVertexData(momentum);
    }

    if constexpr (m > 1 + ST::dim) {
      auto energy = std::make_shared<UMVGF>("E");
      energy->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 1 + ST::dim, 1>(
        "E0", subTriangulation, basis, u, 0));
        energy->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 1 + ST::dim, 1>(
        "E1", subTriangulation, basis, u, 1));
      vtk.addVertexData(energy);
    }

    vtk.write(filename, Dune::VTK::OutputType::appendedraw, 2,Dune::UDGVTKWriteMode::writeVolume);
}

template<class ST, class Vector, class Basis, int m, int finiteElementSize, class ExactSolution>
void writeMultiDomainFunction(const ST& subTriangulation, const Vector& u, const Basis& basis, const std::string& filename,
                              std::integral_constant<int, m>, std::integral_constant<int, finiteElementSize>,
                              const ExactSolution& exactSolution, double t)
{
    using LocalVector = Dune::FieldVector<double, m * finiteElementSize>;

    Dune::RefinedVtkWriter<ST, double> vtk(subTriangulation.gridView(), subTriangulation);

    using UMVGF = UnfittedMultiDomainVTKGridFunction<ST>;
    auto rho = std::make_shared<UMVGF>("rho");
    rho->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 0, 1>(
      "rho0", subTriangulation, basis, u, 0));
    rho->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 0, 1>(
      "rho1", subTriangulation, basis, u, 1));
    vtk.addVertexData(rho);

    if constexpr (m > 1) {
      auto momentum = std::make_shared<UMVGF>("m");
      momentum->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 1, ST::dim>(
        "m0", subTriangulation, basis, u, 0));
      momentum->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 1, ST::dim>(
        "m1", subTriangulation, basis, u, 1));
      vtk.addVertexData(momentum);
    }

    if constexpr (m > 1 + ST::dim) {
      auto energy = std::make_shared<UMVGF>("E");
      energy->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 1 + ST::dim, 1>(
        "E0", subTriangulation, basis, u, 0));
        energy->add(new UnfittedVectorVTKGridFunction<ST, Vector, Basis, LocalVector, 1 + ST::dim, 1>(
        "E1", subTriangulation, basis, u, 1));
      vtk.addVertexData(energy);
    }

    vtk.addVertexData(std::make_shared<ExactSolutionUnfittedVTKGridFunction<ST, m>>(exactSolution, t));

    vtk.write(filename, Dune::VTK::OutputType::appendedraw, 2,Dune::UDGVTKWriteMode::writeVolume);
}

} // end namespace Dune::Hypercut

#endif