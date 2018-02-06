#include <iostream>
#include <fstream>
#include <chrono>
#include <utility>
#include <cmath>
#include "../include/simpleMD.hpp"

using std::vector;
using std::forward_list;
using std::pair;
using std::make_pair;
using std::sqrt;
using std::pow;
using std::exp;
using std::floor;
using uint = unsigned int;

double soft_disk_potential(const vec3& dr) {
    static double rc = pow(2.0, 1.0 / 6.0);
    double dr_mag = dr.mag();
    double potential = 0.0;
    if (dr_mag < rc) {
        potential = 4 * (pow(dr_mag, -12) - pow(dr_mag, -6)) + 1.0;
    }
    return potential;
}

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
                        {-1, 1, 0}, {0, 0, 1}, {1, 0, 1},
                        {1, 1, 1}, {0, 1, 1}, {-1, 1, 1}, {-1, 0, 1},
                        {0, 1, -1}, {-1, 1, -1}, {1, 1, -1}};

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
                          alpha2(exp(-gamma * dt)), 
                          rand_I(0.0, sqrt(1 - pow(alpha2, 2) * temp / mass)),
                          langevin(false), max_v_sum(0.0),
                          n_cells{(uint) floor(box_dim.x / rs), 
                                  (uint) floor(box_dim.y / rs), 
                                  (uint) floor(box_dim.z / rs)},
                          cell_dim{box_dim.x / n_cells[0],
                                   box_dim.y / n_cells[1],
                                   box_dim.z / n_cells[2]} {
    particles.reserve(n_particles);
    nebr_list.reserve(nebr_fac * n_particles);

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

    update_neighbors();
    init_hists();
}

void SimpleMDBox::nvt_step() {
    leapfrog_step();
    time += dt;
}

double SimpleMDBox::velocity_sum() const {
    vec3 sum = {0.0, 0.0, 0.0};
    for (const particle& p : particles) {
        sum += p.v;
    }
    return sum.mag();
}

double SimpleMDBox::kinetic_energy() const {
    double kin = 0.0;
    for (const particle& p : particles) {
        kin += p.v.mag() * p.v.mag();
    }
    return 0.5 * mass * kin;
}

double SimpleMDBox::potential_energy() const {
    double pot = 0.0;
    for (const auto& pair: nebr_list) {
        vec3 dr = pair.p1.r - pair.p2.r;
        wrap_dr(dr);
        pot += soft_disk_potential(dr);
    }

    return pot;
}

void SimpleMDBox::set_langevin(bool enable) {
    langevin = enable;
}

// Updates all particles' accelerations.
void SimpleMDBox::compute_forces() {
    for (auto& p: particles) {
        p.a = {0.0, 0.0, 0.0};
    }

    for (const auto& pair: nebr_list) {
        vec3 dr = pair.p1.r - pair.p2.r;
        wrap_dr(dr);
        vec3 force = soft_disk_force(dr);
        pair.p1.a += force / mass;
        pair.p2.a += -force / mass;
    }
}

void SimpleMDBox::leapfrog_step() {
    static bool update_nebrs = false;

    compute_forces();
    for (particle& p: particles) {
        p.v += 0.5 * p.a * dt;
        p.r += p.v * dt;
    }
    wrap_particles();

    if (update_nebrs) {
        update_nebrs = false;
        max_v_sum = 0.0;
        update_neighbors();
    }

    compute_forces();
    if (langevin) {
        for (particle& p : particles) {
            p.v = alpha2 * (p.v + 0.5 * p.a * dt) + vec3{rand_I(rng), rand_I(rng), rand_I(rng)};
        }       
    } else {
        for (particle& p: particles) {
            p.v += 0.5 * p.a * dt;
        }
    }

    add_max_v();
    if (max_v_sum > 0.5 * nebr_dr) {
        update_nebrs = true;
    }
}

void SimpleMDBox::init_hists() {
    hists.init_hist("v", 0.0, 5.0, 100);
    hists.init_hist("g", 0.8, 1.5, 100);
    hists.set_updater("v", [this](Histogram& hist){
        for (const auto& p : particles) {
            hist.add(p.v.mag());
        }
    });
    hists.set_updater("g", [this](Histogram& hist) {
        for (const auto& pair: nebr_list) {
            vec3 dr = pair.p1.r - pair.p2.r;
            hist.add(dr.mag());
        }
    });
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

void SimpleMDBox::wrap_dr(vec3& dr) const {
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

void SimpleMDBox::add_max_v() {
    double max_vv = 0.0;
    for (const auto& p: particles) {
        const double vv = p.v.square_mag();
        if (vv > max_vv) {
            max_vv = vv;
        }
    }
    max_v_sum += sqrt(max_vv) * dt;
}

void SimpleMDBox::update_neighbors() {
    update_cells();
    nebr_list.clear();
    for (int x = 0; x < n_cells[0]; ++x) {
        for (int y = 0; y < n_cells[1]; ++y) {
            for (int z = 0; z < n_cells[2]; ++z) {
                auto cell = get_cell(x, y, z);
                // Separate case for checking cell against itself
                for (auto p1 = cell.cbegin(); p1 != cell.cend(); ++p1) {
                    auto p2 = p1;
                    ++p2;
                    for (; p2 != cell.cend(); ++p2) {
                        vec3 dr = (*p1)->r - (*p2)->r;
                        wrap_dr(dr);
                        if (dr.mag() < rs) {
                            nebr_list.push_back(particle_pair{**p1, **p2});
                        }
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
                            if (dr.mag() < rs) {
                                nebr_list.push_back(particle_pair{**p1, **p2});
                            }
                        }
                    }
                }
            }
        }
    }
}