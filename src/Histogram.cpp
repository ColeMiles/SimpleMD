#include "../include/Histogram.hpp"
#include <cmath>

Histogram::Histogram(double min_val, double max_val, int num_bins) :
        min_val(min_val), max_val(max_val), histogram(num_bins, {0.0, 0}),
        bin_width( (max_val - min_val) / num_bins ) {
            
    // Label the left edges of the bins
    double prev_edge = min_val;
    for (int bin_num = 0; bin_num < num_bins; ++bin_num) {
        histogram[bin_num].left_edge = prev_edge;
        prev_edge += bin_width;
    }
}

void Histogram::add(double val) {
    int bin_num = (int) std::floor((val - min_val) / bin_width);
    if (bin_num < histogram.size()) {
       ++(histogram[bin_num].count);
    }
}

hist_bin Histogram::get_bin(int bin_num) const {
    return histogram[bin_num];
}

int Histogram::size() const {
    return histogram.size();
}

void Histogram::clear() {
    for (auto& bin: histogram) {
        bin.count = 0;
    }
}

void Histogram::normalize() {
    // Sum up all of the bins
    double sum = 0.0;
    for (const auto& bin : histogram) {
        sum += bin.count;
    }

    // Normalize the area of the histogram to be 1
    const double norm_factor = sum * bin_width;
    // Normalize the histogram
    for (auto& bin: histogram) {
        bin.count /= norm_factor;
    }
}

void Histogram::scale(std::function<double(double)> scale_f) {
    const double half_width = bin_width / 2.0;
    for (auto& bin : histogram) {
        bin.count *= scale_f(bin.left_edge + half_width);
    }
}

std::vector<hist_bin>::iterator Histogram::begin() {
    return histogram.begin();
}

std::vector<hist_bin>::const_iterator Histogram::begin() const {
    return histogram.cbegin();
}

std::vector<hist_bin>::iterator Histogram::end() {
    return histogram.end();
}

std::vector<hist_bin>::const_iterator Histogram::end() const {
    return histogram.cend();
}

