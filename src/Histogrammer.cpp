#include "../include/Histogram.hpp"
#include "../include/Histogrammer.hpp"

void Histogrammer::init_hist(std::string name, double min_val, 
                             double max_val, int num_bins) {
    histograms.emplace(name, Histogram(min_val, max_val, num_bins));
}

void Histogrammer::set_updater(std::string name, 
                               std::function<void(Histogram&)> updater) {
    if (histograms.count(name) == 1) {
        updaters[name] = updater;
    }
}

void Histogrammer::update_hist(std::string name) {
    updaters[name](histograms.at(name));
}

void Histogrammer::clear_hist(std::string name) {
    if (histograms.count(name) == 1) {
        histograms.at(name).clear();
    }
}

Histogram& Histogrammer::operator[](std::string name) {
    return histograms.at(name);
}

Histogram& Histogrammer::get_hist(std::string name) {
    return histograms.at(name);
}

Histogram Histogrammer::get_normalized_hist(std::string name) {
    Histogram copy = histograms.at(name);

    // Sum up all of the bins
    double sum = 0.0;
    for (const auto& bin : copy) {
        sum += bin.count;
    }

    // Normalize the area of the histogram to be 1
    const double norm_factor = sum * copy.bin_width;
    // Normalize the histogram
    for (auto& bin: copy.histogram) {
        bin.count /= norm_factor;
    }

    return copy;
}