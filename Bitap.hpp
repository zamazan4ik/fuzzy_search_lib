#ifndef FUZZYSEARCH_BITAP_HPP
#define FUZZYSEARCH_BITAP_HPP

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>
#include <array>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

template <typename RandomAccessIterator1, typename RandomAccessIterator2>
std::pair<RandomAccessIterator1, RandomAccessIterator1>
bitap_fuzzy_bitwise_search_new(RandomAccessIterator1 corpus_begin, RandomAccessIterator1 corpus_end,
                               RandomAccessIterator2 pattern_begin, RandomAccessIterator2 pattern_end, size_t k)
{
    if (pattern_begin == pattern_end)
    {
        return std::make_pair(corpus_begin, corpus_end);
    }

    typedef typename std::iterator_traits<RandomAccessIterator2>::difference_type PatternType;

    PatternType m = std::distance(pattern_begin, pattern_end);

    //TODO: check sizeof(long) which is platoform-dependent
    assert(m <= 31);

    std::array<unsigned long, std::numeric_limits<char>::max() + 1> pattern_mask;
    pattern_mask.fill(~0);
    for (size_t i = 0; i < m; ++i)
    {
        pattern_mask[pattern_begin[i]] &= ~(1UL << i);
    }

    std::vector<unsigned long> R((k + 1) * sizeof(void*), ~1);
    std::pair<RandomAccessIterator1, RandomAccessIterator1> result(corpus_end, corpus_end);
    for (RandomAccessIterator1 it = corpus_begin; it != corpus_end; ++it)
    {
        /* Update the bit arrays */
        unsigned long old_Rd1 = R[0];

        R[0] |= pattern_mask[*it];
        R[0] <<= 1;

        for (size_t d = 1; d <= k; ++d)
        {
            unsigned long tmp = R[d];
            /* Substitution is all we care about */
            R[d] = (old_Rd1 & (R[d] | pattern_mask[*it])) << 1;
            old_Rd1 = tmp;
        }

        if ((R[k] & (1UL << m)) == 0)
        {
            result.first = (it - m) + 1;
            result.second = result.first + m;
            break;
        }
    }

    return result;
}

template <typename Range1, typename Range2>
std::pair<typename boost::range_iterator<Range1>::type,
          typename boost::range_iterator<Range1>::type>
bitap_fuzzy_bitwise_search_new(Range1 corpus_range, Range2 pattern_range, size_t k)
{
    return bitap_fuzzy_bitwise_search_new(boost::begin(corpus_range), boost::end(corpus_range),
                                          boost::begin(pattern_range), boost::end(pattern_range), k);
}

#endif //FUZZYSEARCH_BITAP_HPP
