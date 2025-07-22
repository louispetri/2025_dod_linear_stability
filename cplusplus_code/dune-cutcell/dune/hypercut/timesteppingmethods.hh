#ifndef DUNE_HYPERCUT_TIMESTEPPING_METHODS_HH
#define DUNE_HYPERCUT_TIMESTEPPING_METHODS_HH

#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>

template<class F, class K_>
class ExplicitEuler
{
public:
    using K = K_;

    ExplicitEuler(F& f) :
        r_(f.zero()),
        f_(f),
        t_(0.0)
    {
    }

    void step(typename F::Domain& x, K deltaT)
    {
        r_ = 0.0;
        f_(t_, x, r_);
        x.axpy(deltaT, r_);
        t_ += deltaT;
    }

private:
    typename F::Domain r_;
    F& f_;
    K t_;
};

template<class F, class K_>
class Heun
{
public:
    using K = K_;

    Heun(F& f) :
        r_(f.zero()),
        y_(f.zero()),
        f_(f),
        t_(0.0)
    {

    }

    void step(typename F::Domain& x, K deltaT)
    {
        r_ = 0.0;
        f_(t_, x, r_);
        y_ = x;
        y_.axpy(deltaT, r_);
        x.axpy(0.5 * deltaT, r_);

        t_ += deltaT;
        r_ = 0.0;
        f_(t_, y_, r_);
        x.axpy(0.5 * deltaT, r_);
    }

private:
    typename F::Domain r_;
    typename F::Domain y_;
    F& f_;
    K t_;
};

template<class F, class K_>
class Shu3
{
public:
    using K = K_;

    Shu3(F& f) :
        r_(f.zero()),
        y_(f.zero()),
        f_(f),
        t_(0.0)
    {

    }

    void step(typename F::Domain& x, K deltaT)
    {
        r_ = 0.0;
        f_(t_, x, r_);
        y_ = x;
        y_.axpy(deltaT, r_);

        r_ = 0.0;
        f_(t_ + deltaT, y_, r_);
        y_ *= 0.25;
        y_.axpy(0.75, x);
        y_.axpy(0.25 * deltaT, r_);

        r_ = 0.0;
        f_(t_ + 0.5 * deltaT, y_, r_);

        x *= (1.0l / 3.0l);
        x.axpy((2.0 / 3.0), y_);
        x.axpy((2.0 / 3.0) * deltaT, r_);

        t_ += deltaT;
    }

private:
    typename F::Domain r_;
    typename F::Domain y_;
    F& f_;
    K t_;
};

template<class F, class K_>
class SSPRK54
{
public:
    using K = K_;

    SSPRK54(F& f) :
        y_{f.zero(), f.zero(), f.zero(), f.zero(), f.zero()},
        r_(f.zero()),
        a_{0.0, 0.0, 0.0, 0.0, 0.0,
           0.39175222700392, 0.0, 0.0, 0.0, 0.0,
           0.21766909633821, 0.36841059262959, 0.0, 0.0, 0.0,
           0.08269208670950, 0.13995850206999, 0.25189177424738, 0.0, 0.0,
           0.06796628370320, 0.11503469844438, 0.20703489864929, 0.54497475021237, 0.0},
        b_{0.14681187618661, 0.24848290924556, 0.10425883036650, 0.27443890091960, 0.22600748319395},
        c_{0.0, 0.39175222700392, 0.58607968896779, 0.47454236302687, 0.93501063100924},
        f_(f),
        t_(0.0)
    {

    }

    void step(typename F::Domain& x, K deltaT)
    {
        for (std::size_t i = 0; i < 5; ++i) {
            y_[i] = 0.0;

            r_ = 0.0;
            r_.axpy(1.0, x);

            for (std::size_t j = 0; j < i; ++j) {
                r_.axpy(deltaT * a_[i][j], y_[j]);
            }

            f_(t_ + deltaT * c_[i], r_, y_[i]);
        }

        for (std::size_t i = 0; i < 5; ++i) {
            x.axpy(deltaT * b_[i], y_[i]);
        }

        t_ += deltaT;
    }

private:
    typename F::Domain y_[5];
    typename F::Domain r_;

    K a_[5][5];
    K b_[5];
    K c_[5];

    F& f_;
    K t_;
};

template<class F, class K_>
class SSPRK104
{
public:
    using K = K_;

    SSPRK104(F& f) :
        y_{f.zero(), f.zero(), f.zero(), f.zero(), f.zero(), f.zero(), f.zero(), f.zero(), f.zero(), f.zero()},
        r_(f.zero()),
        a_{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           1.0 / 6.0, 1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 6.0, 0.0, 0.0, 0.0, 0.0,
           1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 6.0, 1.0 / 6.0, 0.0, 0.0, 0.0,
           1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 0.0, 0.0,
           1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 15.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 0.0},
        b_{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        c_{0.0, 1.0 / 6.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 5.0 / 6.0, 1.0},
        f_(f),
        t_(0.0)
    {

    }

    void step(typename F::Domain& x, K deltaT)
    {
        for (std::size_t i = 0; i < 10; ++i) {
            y_[i] = 0.0;

            r_ = 0.0;
            r_.axpy(1.0, x);

            for (std::size_t j = 0; j < i; ++j) {
                r_.axpy(deltaT * a_[i][j], y_[j]);
            }

            f_(t_ + deltaT * c_[i], r_, y_[i]);
        }

        for (std::size_t i = 0; i < 10; ++i) {
            x.axpy(deltaT * b_[i], y_[i]);
        }

        t_ += deltaT;
    }

private:
    typename F::Domain y_[10];
    typename F::Domain r_;

    K a_[10][10];
    K b_[10];
    K c_[10];

    F& f_;
    K t_;
};

template<class Operator>
std::function<void (typename Operator::Domain&, typename Operator::K)> timesteppingMethod(Operator& op, const std::string& method)
{
  if (method == "heun") {
    auto ts = std::make_shared<Heun<Operator, typename Operator::K>>(op);
    return [ts = std::move(ts)] (typename Operator::Domain& v, typename Operator::K dt) { ts->step(v, dt); };
  } else if (method == "shu3") {
    auto ts = std::make_shared<Shu3<Operator, typename Operator::K>>(op);
    return [ts = std::move(ts)] (typename Operator::Domain& v, typename Operator::K dt) { ts->step(v, dt); };
  } else if (method == "ssprk54") {
    auto ts = std::make_shared<SSPRK54<Operator, typename Operator::K>>(op);
    return [ts = std::move(ts)] (typename Operator::Domain& v, typename Operator::K dt) { ts->step(v, dt); };
  } else if (method == "ssprk104") {
    auto ts = std::make_shared<SSPRK104<Operator, typename Operator::K>>(op);
    return [ts = std::move(ts)] (typename Operator::Domain& v, typename Operator::K dt) { ts->step(v, dt); };
  } else if (method == "expliciteuler") {
    auto ts = std::make_shared<ExplicitEuler<Operator, typename Operator::K>>(op);
    return [ts = std::move(ts)] (typename Operator::Domain& v, typename Operator::K dt) { ts->step(v, dt); };
  } else {
    std::cout << "WARNING: Unknown timestepping method. Reverting to explicit Euler." << std::endl;
    auto ts = std::make_shared<ExplicitEuler<Operator, typename Operator::K>>(op);
    return [ts = std::move(ts)] (typename Operator::Domain& v, typename Operator::K dt) { ts->step(v, dt); };
  }
}

#endif
