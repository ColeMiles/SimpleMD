#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <vector>
#include <utility>
#include <cmath>
#include "../util/vec3.hpp"

using std::vector;
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

// double lj_potential(const vec3& dr) {
//     double dr_mag = dr.mag();
//     return pow(dr_mag, -12) - 2 * pow(dr_mag, -6);
// }

// Book's equation
// rc manually set to 2^(1/6)
double soft_disk_potential(const vec3& dr) {
    static double rc = pow(2.0, 1.0 / 6.0);
    double dr_mag = dr.mag();
    double potential = 0.0;
    if (dr_mag < rc) {
        potential = 4 * (pow(dr_mag, -12) - pow(dr_mag, -6)) + 1.0;
    }
    return potential;
}

// vec3 lj_force(const vec3& dr) {
//     double dr_mag = dr.mag();
//     vec3 force = 12 * (pow(1.0 / dr_mag, 14) - pow(1.0 / dr_mag, 8)) * dr;
//     return force;
// }

// Book's equation
// rc manually set to 2^(1/6)
vec3 soft_disk_force(const vec3& dr) {
    static double rc = pow(2.0, 1.0 / 6.0);
    double dr_mag = dr.mag();
    vec3 force = {0.0, 0.0, 0.0};
    if (dr_mag < rc) {
        force = 48 * (pow(dr_mag, -14) - 0.5 * pow(dr_mag, -8)) * dr;
    }
    return force;
}

// Returns a uniformly-distributed randomly oriented unit vector
vec3 uniform_unit_vec3() {
    static std::mt19937_64 rng(std::chrono::system_clock::now()
                                .time_since_epoch().count());
    static std::normal_distribution<double> uniform;
    vec3 vector = vec3(uniform(rng), uniform(rng), uniform(rng)).normalize();
    return vector;
}

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
    const double rc = pow(2.0, 1.0 / 6.0); // Cut-off distance for LJ potential

    bool langevin;

    std::mt19937_64 rng;
    std::normal_distribution<double> randImpulse;
    
    vec3 compute_net_force(particle& p);
    void wrap_particles();
    void wrap_dr(vec3& dr);
    void compute_forces();
    void leapfrog_step();
    void integrate_positions();
    void integrate_velocities();
    void add_rand_impulses();
};

/* Note to self: As of the moment, the box is entirely in the positive octant,
 *  while in the book it is centered on (0, 0, 0). While both ways should give
 *  the same results, centering on (0, 0, 0) makes a nice symmetry for wrapping
 *  the particle positions both in the position update and in calculating dr
 *  for the force which allows you to use one function for both. Consider 
 *  changing this.
 */

SimpleMDBox::SimpleMDBox(vec3 box_dim, uint n_particles, 
                         double temp, double dt)
                        : box_dim(box_dim), n_particles(n_particles), 
                          temp(temp), dt(dt), time(0.0),
                          volume(box_dim.x * box_dim.y * box_dim.z),
                          rng(std::chrono::system_clock::now()
                                  .time_since_epoch().count()), 
                          randImpulse(0.0, sqrt(2 * gamma * temp * dt / mass)),
                          langevin(false) {
    particles.reserve(n_particles);
    
    // Initialize positions in a simple cubic lattice
    // As of the moment, the entirety of the box is in positive coordinates
    //  (a.k.a. the bottom corner of the box is at (0, 0, 0))
    double cube_vol = volume / n_particles;
    double edge = pow(cube_vol, 1.0 / 3.0);
    int n = 0;

    // Fix to work for general N, and centered at (0.0, 0.0)
    for (double posx = 0.0; posx <= box_dim.x - edge; posx += edge) {
        for (double posy = 0.0; posy <= box_dim.y - edge; posy += edge) {
            for (double posz = 0.0; posz <= box_dim.z - edge; posz += edge) {
                particles.push_back({
                    {posx + edge/2, posy + edge/2, posz + edge/2}, // r
                    {0.0, 0.0, 0.0},                               // v
                    {0.0, 0.0, 0.0}                                // a
                });
                ++n;
            }
        }
    }

    // Velocities initialized as in book, pg. 28
    // (Random directions, magnitude depends on temp, net velocity is zero)
    const double vel_mag = sqrt(3 * (1.0 - 1.0 / n_particles) * temp);
    vec3 net_vel = {0.0, 0.0, 0.0};
    for (particle& p : particles) {
        p.v = vel_mag * uniform_unit_vec3();
        net_vel += p.v;
    }
    const vec3 zero_vel_adjust = -1.0 / n_particles * net_vel;
    for (particle& p : particles) {
        p.v += zero_vel_adjust;
    }

    // Accelerations initialized to zero already when vector was filled
}


// Part of the slow, O(n^2) calc. Will need to be removed when optimizing.
// Also looks a bit sloppy, but since this will be removed it's fine.
vec3 SimpleMDBox::compute_net_force(particle& p) {
    vec3 net_force = {0.0, 0.0, 0.0};
    for (const particle& other: particles) {
        // LJ pair forces
        if (&other != &p) {
            vec3 dr = other.r - p.r;

            // Periodic wrapping
            wrap_dr(dr);
            net_force += soft_disk_force(dr);
        }

        // Damping
        if (langevin) {
            net_force -= gamma * p.v;
        }
    }

    return net_force;
}

// This is different then wrap_dr just because how I set up my
//  coordinates. I should change where the lattice is centered so that this is
//  only one function.
void SimpleMDBox::wrap_particles() {
    for (particle& p : particles) {
        if (p.r.x < 0.0) {
            p.r.x += box_dim.x;
        } else if (p.r.x > box_dim.x) {
            p.r.x -= box_dim.x;
        }
        if (p.r.y < 0.0) {
            p.r.y += box_dim.y;
        } else if (p.r.y > box_dim.y) {
            p.r.y -= box_dim.y;
        }
        if (p.r.z < 0.0) {
            p.r.z += box_dim.z;
        } else if (p.r.z > box_dim.z) {
            p.r.z -= box_dim.z;
        }
    }
}

void SimpleMDBox::wrap_dr(vec3& dr) {
    if (dr.x > box_dim.x / 2) {
        dr.x -= box_dim.x;
    } else if (dr.x < -box_dim.x / 2) {
        dr.x += box_dim.x;
    }
    if (dr.y > box_dim.y / 2) {
        dr.y -= box_dim.y;
    } else if (dr.y < -box_dim.y / 2) {
        dr.y += box_dim.y;
    }
    if (dr.z > box_dim.z / 2) {
        dr.z -= box_dim.z;
    } else if (dr.z < -box_dim.z / 2) {
        dr.z += box_dim.z;
    }
}

// Updates all particles' accelerations.
void SimpleMDBox::compute_forces() {
    for (particle& p: particles) {
        p.a = {0.0, 0.0, 0.0};
        if (langevin) {
            p.a = -gamma * p.v;
        }
    }

    // Loops produce every pair of atoms
    for (auto iter1 = particles.begin(); iter1 != particles.end() - 1; 
                                                                ++iter1) {
        for (auto iter2 = iter1 + 1; iter2 != particles.end(); ++iter2) {
            vec3 dr = iter1->r - iter2->r;
            wrap_dr(dr);
            vec3 force = soft_disk_force(dr);
            iter1->a += force / mass;
            iter2->a += -force / mass;
        }
    }
}

// These two together comprise the VV method
// void SimpleMDBox::integrate_positions() {
//     for (particle& p: particles) {
//         p.r = p.r + p.v * dt + 0.5 * p.a * dt * dt;
//     }
// }

// void SimpleMDBox::integrate_velocities() {
//     for (particle& p: particles) {
//         vec3 new_accel = compute_net_force(p) / mass;
//         p.v += 0.5 * (p.a + new_accel) * dt;
//         p.a = new_accel;
//     }
// }

// Leapfrog method
void SimpleMDBox::leapfrog_step() {
    compute_forces();
    for (particle& p: particles) {
        p.v += 0.5 * p.a * dt;
        p.r += p.v * dt;
    }
    wrap_particles();
    compute_forces();
    for (particle& p: particles) {
        p.v += 0.5 * p.a * dt;
    }
}

void SimpleMDBox::add_rand_impulses() {
    for (particle& p: particles) {
        p.v += randImpulse(rng);
    }
}

// nvt_step for when using VV method
// void SimpleMDBox::nvt_step() {
//     wrap_particles();
//     integrate_positions();
//     integrate_velocities();
//     if (langevin) {
//         add_rand_impulses();
//     }
//     time += dt;
// }

void SimpleMDBox::nvt_step() {
    leapfrog_step();
    if (langevin) {
        add_rand_impulses();
    }
    time += dt;
}

// More of a check for correctness than any physical use
// Adds all the velocity vectors, then takes the magnitude
// (Should always be zero)
double SimpleMDBox::velocity_sum() {
    vec3 sum = {0.0, 0.0, 0.0};
    for (particle&p : particles) {
        sum += p.v;
    }
    return sum.mag();
}

double SimpleMDBox::kinetic_energy() {
    double kin = 0.0;
    for (particle& p : particles) {
        kin += 0.5 * mass * p.v.mag() * p.v.mag();
    }
    return kin;
}

double SimpleMDBox::potential_energy() {
    double pot = 0.0;
    for (auto iter1 = particles.cbegin(); iter1 != particles.cend() - 1; 
                                                                 ++iter1) {
        for (auto iter2 = iter1 + 1; iter2 != particles.cend(); ++iter2) {
            vec3 dr = iter1->r - iter2->r;
            
            // Periodic wrapping
            wrap_dr(dr);

            pot += soft_disk_potential(dr);
        }
    }
    return pot;
}

double SimpleMDBox::total_energy() {
    return kinetic_energy() + potential_energy();
}

// Returns a vector of tuples, where each tuple is (bin_left_edge, rel_freq)
vector<pair<double,double>> SimpleMDBox::compute_velocity_hist(double min_val, 
                                                                double max_val,
                                                                int num_bins) {
    double bin_width = (max_val - min_val) / num_bins;
    vector<pair<double,double>> histogram(num_bins, {0.0, 0.0});

    // Label the left edges of the bins
    double prev_edge = 0.0;
    for (int bin = 1; bin < num_bins; ++bin) {
        prev_edge += bin_width;
        histogram[bin].first = prev_edge;
    }

    // Make a tally of the number of particles in each bin
    for (particle& p : particles) {
        double vel = p.v.mag();
        long bin = (long) floor(vel / bin_width);
        if (bin < num_bins) {
           ++(histogram[bin].second);
        }
    }

    // Normalize the histogram
    for (int bin = 0; bin < num_bins; ++bin) {
        histogram[bin].second /= n_particles;
    }

    return histogram;
}

void SimpleMDBox::set_langevin(bool enable) {
    langevin = enable;
}