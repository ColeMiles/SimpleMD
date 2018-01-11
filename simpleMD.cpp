#include <iostream>
#include <fstream>
#include <chrono>
#include <utility>
#include <cmath>
#include "simpleMD.hpp"

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

// Note a difference from the book: I am not including {0, 0, 0} as an offset,
//  that case is handled separately
const vec3 SimpleMDBox::OFFSETS[] = {{1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                        {-1, 1, 0}, {-1, 0, 0}, {0, 0, 1}, {1, 0, 1},
                        {1, 1, 1}, {0, 1, 1}, {-1, 1, 1}, {-1, 0, 1},
                        {0, 1, -1}, {-1, 1, -1}};

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
                          rand_impulse(0.0, sqrt(2 * gamma * temp * dt / mass)),
                          langevin(false),
                          n_cells{(uint) floor(box_dim.x / rc), 
                                  (uint) floor(box_dim.y / rc), 
                                  (uint) floor(box_dim.z / rc)},
                          cell_dim{box_dim.x / n_cells[0],
                                   box_dim.y / n_cells[1],
                                   box_dim.z / n_cells[2]} {
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

    // Initialize the cell list. The number of cells in each dimension is found
    //  as the number of cells to make the edge length just barely more than rc
    //  (this is done in the initialization list at the top)
    update_cells();
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

inline forward_list<particle*>& SimpleMDBox::get_cell(int x, int y, int z) {
    static const uint offset_y = n_cells[1];
    static const uint offset_z = n_cells[1] * n_cells[2];
    return cell_list[x + y * offset_y + z * offset_z];
}

inline void SimpleMDBox::wrap_cell(int& x, int& y, int& z) {
    if (x >= (int) n_cells[0]) {
        x -= n_cells[0];
    } else if (x < 0) {
        x += n_cells[0];
    }
    if (y >= (int) n_cells[1]) {
        y -= n_cells[1];
    } else if (y < 0) {
        y += n_cells[1];
    }
    if (z >= (int) n_cells[2]) {
        z -= n_cells[2];
    } else if (z < 0) {
        z += n_cells[2];
    }
}

/* Updates all of the lists of particles in each cell.
 * Always call this after wrapping all of the particles, never before
 */
void SimpleMDBox::update_cells() {
    // These offsets are used because we need our 1D vector to act like a
    //   3D matrix, so these offsets are used to calculate the index 
    //   using 3 coordinates
    uint offset_y = n_cells[1];
    uint offset_z = n_cells[1] * n_cells[2];
    // Need a new vector of empty lists
    cell_list = vector<forward_list<particle*>>(n_cells[0] * n_cells[1]
                                                           * n_cells[2]);
    for (particle & p : particles) {
        uint cell_x = floor(p.r.x / cell_dim.x);
        uint cell_y = floor(p.r.y / cell_dim.y);
        uint cell_z = floor(p.r.z / cell_dim.z);
        auto& cell = get_cell(cell_x, cell_y, cell_z);
        cell.emplace_front(&p);
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

    // // Loops produce every pair of atoms
    // for (auto iter1 = particles.begin(); iter1 != particles.end() - 1; 
    //                                                             ++iter1) {
    //     for (auto iter2 = iter1 + 1; iter2 != particles.end(); ++iter2) {
    //         vec3 dr = iter1->r - iter2->r;
    //         wrap_dr(dr);
    //         vec3 force = soft_disk_force(dr);
    //         iter1->a += force / mass;
    //         iter2->a += -force / mass;
    //     }
    // }

    // Loop over all cells, then over all the neighboring cells specified by
    //  OFFSETS
    for (int x = 0; x < n_cells[0]; ++x) {
        for (int y = 0; y < n_cells[1]; ++y) {
            for (int z = 0; z < n_cells[2]; ++z) {
                auto cell = get_cell(x, y, z);
                // Separate case for checking cell against itself
                for (auto p1 = cell.cbegin(); p1 != cell.cend(); ++p1) {
                    auto copy = p1;
                    for (auto p2 = ++copy; p2 != cell.cend(); ++p2) {
                        vec3 dr = (*p1)->r - (*p2)->r;
                        wrap_dr(dr);
                        vec3 force = soft_disk_force(dr);
                        (*p1)->a += force / mass;
                        (*p2)->a += -force / mass;
                    }
                }

                for (const vec3& offset: OFFSETS) {
                    int nx = x + offset.x;
                    int ny = y + offset.y;
                    int nz = z + offset.z;
                    wrap_cell(nx, ny, nz);
                    auto neighbor = get_cell(nx, ny, nz);
                    // Loop over all pairs of particles in these two cells, do
                    //  as before
                    for (auto p1 = cell.cbegin(); p1 != cell.cend(); ++p1) {
                        for (auto p2 = neighbor.cbegin(); p2 != neighbor.cend(); 
                                                                         ++p2) {
                            vec3 dr = (*p1)->r - (*p2)->r;
                            wrap_dr(dr);
                            vec3 force = soft_disk_force(dr);
                            (*p1)->a += force / mass;
                            (*p2)->a += -force / mass;
                        }
                    }
                }
            }
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
    update_cells();
    compute_forces();
    for (particle& p: particles) {
        p.v += 0.5 * p.a * dt;
    }
}

void SimpleMDBox::add_rand_impulses() {
    for (particle& p: particles) {
        p.v += rand_impulse(rng);
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