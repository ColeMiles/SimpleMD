#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>

struct hist_bin {
        double left_edge;
        double count;
};

// Represents a single histogram of some data
class Histogram {
friend class Histogrammer;

public:
    Histogram(double min_val, double max_val, int num_bins);
    void add(double val);
    hist_bin get_bin(int bin_num) const;
    int size() const;
    void clear();

    // Makes the class work with range-based for loops
    std::vector<hist_bin>::iterator begin();
    std::vector<hist_bin>::const_iterator begin() const;
    std::vector<hist_bin>::iterator end();
    std::vector<hist_bin>::const_iterator end() const;
private:
    double bin_width, min_val, max_val;
    std::vector<hist_bin> histogram;
};

#endif