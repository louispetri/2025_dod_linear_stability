#ifndef DUNE_SUBTRIANGULATION_EXAMPLEDOMAINS_HH
#define DUNE_SUBTRIANGULATION_EXAMPLEDOMAINS_HH

#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/subtriangulation/simpletpmctriangulation/domainconfiguration.hh>

namespace Dune::SubTriangulation {

template<class K_>
class Ramp
{
public:
    static const int dim = 2;
    using K = K_;
    using Coordinate = Dune::FieldVector<K, dim>;

    using Grid = Dune::YaspGrid<dim>;
    using GridView = typename Grid::LevelGridView;

    Ramp(int N, K cutAngle, K cutOffset, K backgroundSize = 1.0)
        : grid_({backgroundSize, backgroundSize}, {N, N})
        , gridView_(grid_.levelGridView(0))
        , levelGridView_(grid_.levelGridView(grid_.maxLevel()))
    {
        K h =  1.0l / N;
        Dune::FieldVector<K, dim> normal { std::cos(cutAngle), std::sin(cutAngle) };

        auto levelset = [normal, cutOffset](Dune::FieldVector<K, dim> x) -> K
        {
            x[0] -= cutOffset;
            return -normal * x;
        };

        domainConfiguration_.addInterface(Interface<GridView>(0, Dune::Functions::makeGridViewFunction(levelset, levelGridView_)));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(0, {"i"}));
    }

    const GridView& gridView() const
    {
        return gridView_;
    }

    const GridView& levelGridView() const
    {
        return levelGridView_;
    }

    const DomainConfiguration<GridView, GridView>& domainConfiguration() const
    {
        return domainConfiguration_;
    }

    K backgroundEdgeLength() const
    {
        return gridView_.template begin<1>()->geometry().volume();
    }

private:
    Grid grid_;
    GridView gridView_;
    GridView  levelGridView_;
    DomainConfiguration<GridView, GridView> domainConfiguration_;
};

template<class K_>
class UnitSquareWithStraightCut
{
public:
    static const int dim = 2;
    using K = K_;
    using Coordinate = Dune::FieldVector<K, dim>;

    using Grid = Dune::YaspGrid<dim>;
    using GridView = typename Grid::LevelGridView;

    UnitSquareWithStraightCut(int N, K cutAngle, K cutOffset, K backgroundSize = 1.0l)
        : grid_({backgroundSize, backgroundSize}, {N , N})
        , gridView_(grid_.levelGridView(0))
        , levelGridView_(grid_.levelGridView(grid_.maxLevel()))
    {
        K h =  1.0l / N;
        K angle = (cutAngle / 180.0l) * M_PI + 0.5 * M_PI;
        Dune::FieldVector<K, dim> normal { std::cos(angle), std::sin(angle) };

        auto levelset = [normal, cutOffset](Dune::FieldVector<K, dim> x) -> K 
        {
            x[0] -= cutOffset;
            return -normal * x;
        };

        domainConfiguration_.addInterface(Interface<GridView>(0, Dune::Functions::makeGridViewFunction(levelset, levelGridView_)));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(0, {"i"}));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(1, {"e"}));
    }

    const GridView& gridView() const
    {
        return gridView_;
    }

    const GridView& levelGridView() const
    {
        return levelGridView_;
    }

    const DomainConfiguration<GridView, GridView>& domainConfiguration() const
    {
        return domainConfiguration_;
    }

    K backgroundEdgeLength() const
    {
        return gridView_.template begin<1>()->geometry().volume();
    }

private:
    Grid grid_;
    GridView gridView_;
    GridView  levelGridView_;
    DomainConfiguration<GridView, GridView> domainConfiguration_;
};

template<class K_>
class Disk
{
public:
    static const int dim = 2;
    using K = K_;
    using Coordinate = Dune::FieldVector<K, dim>;

    using Grid = Dune::YaspGrid<dim>;
    using GridView = typename Grid::LevelGridView;

    Disk(int N, K radius, K centerOffsetX)
        : grid_({2 * radius + 0.1, 2 * radius + 0.1}, {N, N})
        , gridView_(grid_.levelGridView(0))
        , levelGridView_(grid_.levelGridView(grid_.maxLevel()))
        , center_(radius + 0.05)
        , radius_(radius)
    {
        center_[0] += centerOffsetX;

        auto levelset = [this](Dune::FieldVector<K, dim> x) -> K
        {
            x -= center_;
            return x.two_norm() - radius_;
        };

        domainConfiguration_.addInterface(Interface<GridView>(0, Dune::Functions::makeGridViewFunction(levelset, levelGridView_)));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(0, {"i"}));
    }

    const GridView& gridView() const
    {
        return gridView_;
    }

    const GridView& levelGridView() const
    {
        return levelGridView_;
    }

    const DomainConfiguration<GridView, GridView>& domainConfiguration() const
    {
        return domainConfiguration_;
    }

    Coordinate center() const
    {
        return center_;
    }

    K radius() const
    {
        return radius_;
    }

    K backgroundEdgeLength() const
    {
        return gridView_.template begin<1>()->geometry().volume();
    }

private:
    Grid grid_;
    GridView gridView_;
    GridView  levelGridView_;
    DomainConfiguration<GridView, GridView> domainConfiguration_;
    Coordinate center_;
    K radius_;
};

template<class K_>
class RotatedUnitSquare
{
public:
    static const int dim = 2;
    using K = K_;
    using Coordinate = Dune::FieldVector<K, dim>;

    using Grid = Dune::YaspGrid<dim>;
    using GridView = typename Grid::LevelGridView;

    using ReferenceGeometry = Dune::CachedMultiLinearGeometry<K, dim, dim>;

    RotatedUnitSquare(int N, K cutAngle)
        : grid_(Coordinate({std::cos(cutAngle) + std::sin(cutAngle), std::cos(cutAngle) + std::sin(cutAngle)}), {N, N})
        , gridView_(grid_.levelGridView(0))
        , levelGridView_(grid_.levelGridView(grid_.maxLevel()))
    {
        x_offset_ = std::sin(cutAngle);
        y_offset_ = std::cos(cutAngle);

        K angle = cutAngle;
        Dune::FieldVector<K, dim> normal { std::cos(angle), std::sin(angle) };

        Coordinate lowerPoint;
        lowerPoint[0] = x_offset_;
        lowerPoint[1] = 0.0;

        auto levelset = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= lowerPoint;
            return -normal * x;
        };

        angle += 0.5 * M_PI;
        Dune::FieldVector<K, dim> normal2 { std::cos(angle), std::sin(angle) };

        auto levelset2 = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= lowerPoint;
            return -normal2 * x;
        };

        angle += 0.5 * M_PI;
        Dune::FieldVector<K, dim> normal3 { std::cos(angle), std::sin(angle) };

        lowerPoint[0] = x_offset_ + (1.0 / y_offset_);

        auto levelset3 = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= lowerPoint;
            return -normal3 * x;
        };

        angle += 0.5 * M_PI;
        Dune::FieldVector<K, dim> normal4 { std::cos(angle), std::sin(angle) };

        lowerPoint[0] = -1.0 / x_offset_ + x_offset_;

        auto levelset4 = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= lowerPoint;
            return -normal4 * x;
        };

        domainConfiguration_.addInterface(Interface<GridView>(0, Dune::Functions::makeGridViewFunction(levelset, levelGridView_)));
        domainConfiguration_.addInterface(Interface<GridView>(1, Dune::Functions::makeGridViewFunction(levelset2, levelGridView_)));
        domainConfiguration_.addInterface(Interface<GridView>(2, Dune::Functions::makeGridViewFunction(levelset3, levelGridView_)));
        domainConfiguration_.addInterface(Interface<GridView>(3, Dune::Functions::makeGridViewFunction(levelset4, levelGridView_)));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(0, {"iiii"}));
    }

    const GridView& gridView() const
    {
        return gridView_;
    }

    const GridView& levelGridView() const
    {
        return levelGridView_;
    }

    const DomainConfiguration<GridView, GridView>& domainConfiguration() const
    {
        return domainConfiguration_;
    }

    K backgroundEdgeLength() const
    {
        return gridView_.template begin<1>()->geometry().volume();
    }

    K xOffset() const
    {
        return x_offset_;
    }

    K yOffset() const
    {
        return y_offset_;
    }

    ReferenceGeometry referenceGeometry() const
    {
        std::vector<Coordinate> points;
        points.resize(4, Coordinate(0.0));

        points[0][0] = x_offset_;

        points[1][0] = x_offset_ + y_offset_;
        points[1][1] = x_offset_;

        points[2][1] = y_offset_;

        points[3][0] = y_offset_;
        points[3][1] = x_offset_ + y_offset_;

        return ReferenceGeometry(Dune::GeometryTypes::cube(dim), points);
    }

private:
    Grid grid_;
    GridView gridView_;
    GridView  levelGridView_;
    DomainConfiguration<GridView, GridView> domainConfiguration_;
    K x_offset_;
    K y_offset_;
};

template<class K_>
class RotatedOpenChannel
{
public:
    static const int dim = 2;
    using K = K_;
    using Coordinate = Dune::FieldVector<K, dim>;

    using Grid = Dune::YaspGrid<dim, EquidistantOffsetCoordinates<K, dim>>;
    using GridView = typename Grid::LevelGridView;

    RotatedOpenChannel(int N, K length, K cutAngle, K offset)
        : grid_(Coordinate({0.0, 0.0}), Coordinate({length, length}), {N, N})
        , gridView_(grid_.levelGridView(0))
        , levelGridView_(grid_.levelGridView(grid_.maxLevel()))
    {
        K h =  1.0l / N;
        K angle = (cutAngle / 180.0l) * M_PI + 0.5 * M_PI;
        Dune::FieldVector<K, dim> normal { std::cos(angle), std::sin(angle) };

        Coordinate lowerPoint;
        lowerPoint[0] = offset;
        lowerPoint[1] = 0.0;

        auto levelset = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= lowerPoint;
            return -normal * x;
        };

        lowerPoint[0] = 0.0;
        lowerPoint[1] = offset;

        auto levelset2 = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= lowerPoint;
            return normal * x;
        };

        domainConfiguration_.addInterface(Interface<GridView>(0, Dune::Functions::makeGridViewFunction(levelset, levelGridView_)));
        domainConfiguration_.addInterface(Interface<GridView>(1, Dune::Functions::makeGridViewFunction(levelset2, levelGridView_)));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(0, {"ii"}));
    }

    const GridView& gridView() const
    {
        return gridView_;
    }

    const GridView& levelGridView() const
    {
        return levelGridView_;
    }

    const DomainConfiguration<GridView, GridView>& domainConfiguration() const
    {
        return domainConfiguration_;
    }

    K backgroundEdgeLength() const
    {
        return gridView_.template begin<1>()->geometry().volume();
    }

private:
    Grid grid_;
    GridView gridView_;
    GridView  levelGridView_;
    DomainConfiguration<GridView, GridView> domainConfiguration_;
};

template<class K_>
class WitchOfAgnesi
{
public:
    static const int dim = 2;
    using K = K_;
    using Coordinate = Dune::FieldVector<K, dim>;

    using Grid = Dune::YaspGrid<dim, EquidistantOffsetCoordinates<K, dim>>;
    using GridView = typename Grid::LevelGridView;

    WitchOfAgnesi(int xCells, int yCells, Coordinate& lowerLeft, Coordinate& upperRight, K a, K h0)
        : grid_(lowerLeft, upperRight, {xCells, yCells})
        , gridView_(grid_.levelGridView(0))
        , levelGridView_(grid_.levelGridView(grid_.maxLevel()))

    {
        auto levelset = [=](Dune::FieldVector<K, dim> x) -> K
        {
            return (h0 * a * a) / (x[0] * x[0] + a * a) - x[1];
        };

        domainConfiguration_.addInterface(Interface<GridView>(0, Dune::Functions::makeGridViewFunction(levelset, levelGridView_)));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(0, {"i"}));
    }

    const GridView& gridView() const
    {
        return gridView_;
    }

    const GridView& levelGridView() const
    {
        return levelGridView_;
    }

    const DomainConfiguration<GridView, GridView>& domainConfiguration() const
    {
        return domainConfiguration_;
    }

    K backgroundEdgeLength() const
    {
        return gridView_.template begin<1>()->geometry().volume();
    }

private:
    Grid grid_;
    GridView gridView_;
    GridView  levelGridView_;
    DomainConfiguration<GridView, GridView> domainConfiguration_;
};

template<class K_>
class PeriodicRectangleWithStraightCut
{
public:
    static const int dim = 2;
    using K = K_;
    using Coordinate = Dune::FieldVector<K, dim>;

    using Grid = Dune::YaspGrid<dim>;
    using GridView = typename Grid::LevelGridView;

    PeriodicRectangleWithStraightCut(int N, int M, K cutAngle1, K cutOffset1, K cutAngle2, K cutOffset2, K width, K height)
        : grid_({width, height}, {N , M}, std::bitset<dim>("01"), 1, MPI_COMM_WORLD)
        , gridView_(grid_.levelGridView(0))
        , levelGridView_(grid_.levelGridView(grid_.maxLevel()))
    {
        K angle1 = (cutAngle1 / 180.0l) * M_PI + 0.5 * M_PI;
        Dune::FieldVector<K, dim> normal1 { std::cos(angle1), std::sin(angle1) };

        auto levelset = [normal1, cutOffset1](Dune::FieldVector<K, dim> x) -> K
        {
            x[0] -= cutOffset1;
            return -normal1 * x;
        };

        K angle2 = (cutAngle2 / 180.0l) * M_PI + 0.5 * M_PI;
        Dune::FieldVector<K, dim> normal2 { std::cos(angle2), std::sin(angle2) };

        auto levelset2 = [normal2, cutOffset2](Dune::FieldVector<K, dim> x) -> K
        {
            x[0] -= cutOffset2;
            return normal2 * x;
        };

        domainConfiguration_.addInterface(Interface<GridView>(0, Dune::Functions::makeGridViewFunction(levelset, levelGridView_)));
        domainConfiguration_.addInterface(Interface<GridView>(1, Dune::Functions::makeGridViewFunction(levelset2, levelGridView_)));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(0, {"ie", "ei"}));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(1, {"ee"}));
    }

    const GridView& gridView() const
    {
        return gridView_;
    }

    const GridView& levelGridView() const
    {
        return levelGridView_;
    }

    const DomainConfiguration<GridView, GridView>& domainConfiguration() const
    {
        return domainConfiguration_;
    }

    K backgroundEdgeLength() const
    {
        return gridView_.template begin<1>()->geometry().volume();
    }

private:
    Grid grid_;
    GridView gridView_;
    GridView  levelGridView_;
    DomainConfiguration<GridView, GridView> domainConfiguration_;
};

template<class K_>
class RotatedPeriodicChannel
{
public:
    static const int dim = 2;
    using K = K_;
    using Coordinate = Dune::FieldVector<K, dim>;

    using Grid = Dune::YaspGrid<dim>;
    using GridView = typename Grid::LevelGridView;

    RotatedPeriodicChannel(int N, K x1, K x2, K y1, K y2)
        : grid_(Coordinate({1.0, 1.0}), {N, N}, std::bitset<dim>("11"), 1, MPI_COMM_WORLD)
        , gridView_(grid_.levelGridView(0))
        , levelGridView_(grid_.levelGridView(grid_.maxLevel()))
    {
        Coordinate normal1{-y1, 1.0 -x2};
        Coordinate normal2{-y2, 1.0 -x1};
        Coordinate normal3{y1 - 1.0, x2};
        Coordinate normal4{y2 - 1.0, x1};

        Coordinate offset1{x2, 0.0};
        Coordinate offset2{x1, 0.0};
        Coordinate offset3{(-y1 * x2) / (1 - y1), 0.0};
        Coordinate offset4{(-y2 * x1) / (1 - y2), 0.0};

        auto levelset = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= offset1;
            return -normal1 * x;
        };

        auto levelset2 = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= offset2;
            return -normal2 * x;
        };

        auto levelset3 = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= offset3;
            return -normal3 * x;
        };

        auto levelset4 = [=](Dune::FieldVector<K, dim> x) -> K
        {
            x -= offset4;
            return -normal4 * x;
        };

        domainConfiguration_.addInterface(Interface<GridView>(0, Dune::Functions::makeGridViewFunction(levelset, levelGridView_)));
        domainConfiguration_.addInterface(Interface<GridView>(1, Dune::Functions::makeGridViewFunction(levelset2, levelGridView_)));
        domainConfiguration_.addInterface(Interface<GridView>(0, Dune::Functions::makeGridViewFunction(levelset3, levelGridView_)));
        domainConfiguration_.addInterface(Interface<GridView>(1, Dune::Functions::makeGridViewFunction(levelset4, levelGridView_)));
        domainConfiguration_.addDomain(Dune::SubTriangulation::Domain(0, {"ieee", "iiie"}));
    }

    const GridView& gridView() const
    {
        return gridView_;
    }

    const GridView& levelGridView() const
    {
        return levelGridView_;
    }

    const DomainConfiguration<GridView, GridView>& domainConfiguration() const
    {
        return domainConfiguration_;
    }

    K backgroundEdgeLength() const
    {
        return gridView_.template begin<1>()->geometry().volume();
    }

private:
    Grid grid_;
    GridView gridView_;
    GridView  levelGridView_;
    DomainConfiguration<GridView, GridView> domainConfiguration_;
};

}

#endif