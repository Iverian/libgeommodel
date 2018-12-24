#ifndef GEOM_MODEL_SRC_GEOM_COX_DE_BOOR_H_
#define GEOM_MODEL_SRC_GEOM_COX_DE_BOOR_H_

#include <geom_model/point.h>

#include <vector>

size_t get_interval(double u, const std::vector<double>& knots);

std::pair<Point, double> cox_de_boor(double u, size_t p, size_t order,
                                     const std::vector<double>& knots,
                                     const std::vector<Point>& points,
                                     const std::vector<double>& weights);

template <class T>
struct VectorView {
    using value_type = typename std::vector<T>::value_type;
    using size_type = typename std::vector<T>::size_type;
    using const_reference = typename std::vector<T>::const_reference;
    using const_pointer = typename std::vector<T>::const_pointer;

    const_reference operator[](size_type i) const
    {
        return data_[i];
    }

    const_reference front() const
    {
        return data_[0];
    }

    const_reference back() const
    {
        return data_[size_ - 1];
    }

    size_type size() const
    {
        return size_;
    }

    VectorView()
        : data_(nullptr)
        , size_(0)
    {
    }

    VectorView(const std::vector<T>& source)
        : data_(source.data())
        , size_(source.size())
    {
    }

    const_pointer data_;
    size_type size_;
};

struct CoxDeBoor {
    // using ResultType = std::pair<Point, double>;
    struct ResultType {
        Point r;
        double w;
    };

    struct CoxDeBoorProxy {
        friend struct CoxDeBoor;
        ResultType get(size_t k) const;
        std::vector<ResultType> range(size_t k) const;

        Point f() const;
        Vec df() const;
        Vec df2() const;

        double u;
        size_t p;

    private:
        CoxDeBoorProxy(double pu, size_t pp, const CoxDeBoor& parent);
        const CoxDeBoor& parent_;
    };

    CoxDeBoor();
    CoxDeBoor(size_t p_order, const std::vector<double>& p_knots,
              const std::vector<Point>& p_points,
              const std::vector<double>& p_weights);

    CoxDeBoorProxy get_proxy(double u) const;

    size_t order;
    VectorView<double> knots;
    VectorView<Point> points;
    VectorView<double> weights;

protected:
    double pad(double u) const;
    size_t get_interval(double u) const;
    ResultType eval(double u, size_t p, std::vector<double>& w,
                    std::vector<Point>& wr) const;
    double get_w(size_t i, size_t k) const;
    Point get_wr(size_t i, size_t k) const;

private:
};

#endif // GEOM_MODEL_SRC_GEOM_COX_DE_BOOR_H_