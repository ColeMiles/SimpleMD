#ifndef SIMPLEMD_HPP
#define SIMPLEMD_HPP

#include <random>
#include <vector>
#include <forward_list>
#include "../util/vec3.hpp"

using std::vector;
using std::forward_list;
using std::pair;
using std::make_pair;
using std::sqrt;
using std::pow;
using std::floor;
using uint = unsigned int;

struct particle {
    vec3 r, v, a;
};

// Unused at the moment
struct particle_pair {
    particle& a;
    particle& b;
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
    
    const uint n_particles;
    double time;

    void nvt_step();
    double velocity_sum();
    double kinetic_energy();
    double potential_energy();
    double total_energy();

    vector<pair<double,double>> compute_velocity_hist(double min_val,
                                                      double max_val,
                                                      int num_bins);

    void set_langevin(bool enable);
private:
    vector<particle> particles;

    const vec3 box_dim;
    const double volume;

    const double temp;
    const double mass = 1.0;
    const double gamma = 0.5;
    const double dt;

    // Cut-off distance for the soft-sphere potential
    constexpr static double rc = pow(2.0, 1.0 / 6.0);

    bool langevin;

    std::mt19937_64 rng;
    std::normal_distribution<double> rand_impulse;
    
    static const vec3 OFFSETS[];
    // Want this to be a vec3, but that's doubles (make it a template?)
    // This is the number of cells in each direction (x, y, z)
    const uint n_cells[3];
    // This is the dimensions of any single cell
    const vec3 cell_dim;
    vector<forward_list<particle*>> cell_list;

    vec3 compute_net_force(particle& p);
    void wrap_particles();
    void wrap_dr(vec3& dr);
    forward_list<particle*>& get_cell(int x, int y, int z);
    // Not pretty - again would be fixed by templating vec3
    void wrap_cell(int& x, int& y, int& z);
    void update_cells();
    void compute_forces();
    void leapfrog_step();
    void integrate_positions();
    void integrate_velocities();
    void add_rand_impulses();
};

#endif