#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include <map>
#include <functional>
#include "Histogram.hpp"

// Manages multiple histograms, and updating them with chunks of data
class Histogrammer {
public:
    // Sets up a histogram in the list of tracked histograms.
    // If a histogram of this name already exists, it will be overwritten.
    void init_hist(std::string name, double min_val, double max_val,
                                                          int num_bins);

    // Sets a function that will be called whenever update_hist is called
    //  with that name. If one already exists with that name, overwrites
    //  it. If the matching hist doesn't exist yet, does nothing.
    void set_updater(std::string name, std::function<void(Histogram&)> updater);
    
    // Calls the function set up to be the updater
    void update_hist(std::string name);

    // Empties a histogram
    void clear_hist(std::string name);

    // Getters to look at the histogram of the given name 
    Histogram& operator[](std::string name);   //<-------- equivalent
    Histogram& get_hist(std::string name);     //<-----|

private:
    std::map<std::string, Histogram> histograms;
    std::map<std::string, std::function<void(Histogram&)>> updaters;
};

#endif