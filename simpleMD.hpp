#ifndef SIMPLEMD_HPP
#define SIMPLEMD_HPP

#include <random>
#include <vector>
#include <forward_list>
#include "util/vec3.hpp"

struct particle {
    vec3 r, v, a;
};

// Unused at the moment
struct particle_pair {
    particle& p1;
    particle& p2;
};

/* NOTES:
 * As of the moment, only works if you
 *  1) Make all dimensions of the box the same (a cube)
 *  2) Choose n_particles to be a perfect cube
 *
 * Which indicates that there's work to be done on the lattice initialization
 */ 

class SimpleMDBox {
public:
    SimpleMDBox(vec3 box_dim, uint N, double temp, double dt);
    
    const unsigned int n_particles;
    double time;

    void nvt_step();
    double velocity_sum();
    double kinetic_energy();
    double potential_energy();

    std::vector<std::pair<double,double>> compute_velocity_hist(double min_val,
                                                                double max_val,
                                                                int num_bins);

    void set_langevin(bool enable);
private:
    std::vector<particle> particles;

    const vec3 box_dim;
    const double volume;

    const double temp;
    const double mass = 1.0;
    const double gamma = 0.5;
    const double alpha2;
    const double dt;

    //=======Cut-off distance for the soft-sphere potential==========
    constexpr static double rc = std::pow(2.0, 1.0 / 6.0);

    //=======Integration functions============
    void compute_forces();
    void leapfrog_step();

    //=======Periodic boundary condition functions==========
    void wrap_particles();
    void wrap_dr(vec3& dr);

    //=======Used in Langevin integration=============
    bool langevin;
    std::mt19937_64 rng;
    std::normal_distribution<double> rand_I;
    
    //=======Neighbor List management===========
    static const vec3 OFFSETS[];
    // The extra length added to rc to get the range to be a "neighbor"
    constexpr static double nebr_dr = 0.3;
    constexpr static double rs = rc + nebr_dr;
    // A place to keep a running sum of the max vel
    double max_v_sum;
    void add_max_v();
    // A guess at the maximum neighbors per particle we'll have
    const unsigned int nebr_fac = 8;
    std::vector<particle_pair> nebr_list;

    // Want this to be a vec3, but that's doubles (make it a template?)
    // This is the number of cells in each direction (x, y, z)
    const unsigned int n_cells[3];
    // This is the dimensions of any single cell
    const vec3 cell_dim;
    std::vector<std::forward_list<particle*>> cell_list;

    std::forward_list<particle*>& get_cell(int x, int y, int z);
    // Not pretty - again would be fixed by templating vec3
    void wrap_cell(int& x, int& y, int& z);
    void update_cells();
    void update_neighbors();
};

#endif
