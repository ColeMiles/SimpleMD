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