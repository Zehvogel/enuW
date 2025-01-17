#include <concepts>

// size_t for N is horrible overkill but :shrug:
template<std::convertible_to<double> T, std::size_t N>
struct Covariance {
    std::size_t m_n = N;
    std::vector<double> m_means{N, 1.0};
    double m_sum_w;

    Covariance()
    {
        // TODO: can probably also live with 1 if I don't make any stupid mistakes in the logic...
        std::static_assert(n >= 2, "Covariance needs at least two columns")
    }

    Covariance(Covariance&) = default;

    void Merge(const std::vector<Covariance>& other)
    {
        m_means = m_sum_w * m_means + other.m_sum_w * other * m_means;
        m_sum_w += other.m_sum_w;
        m_means /= m_sum_w;
    }

    template<T...>
    void Fill(T... args, T last)
    {
        // n+1 columns means last column will be treated as weight
        std::static_assert(sizeof...(args) == m_n || sizeof...(args) == m_n - 1, "Tried to fill with more columns than initialised");
        if (sizeof...(args) == m_n)
            FillWeighted({args...}, last);
        else
            FillWeighted({args..., last}, 1.0);
    }

    private:
    void FillWeighted(std::vector<double> args, double weight)
    {
        m_sum_w += weight;
        for (int i = 0; i <  m_n; i++) {
            m_means[i] += (weight / m_sum_w) * (args[i] - m_means[i]);
        }

    }

};