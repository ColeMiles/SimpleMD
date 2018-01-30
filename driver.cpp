#include "simpleMD.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

using std::cout;
using std::ofstream;
using std::setw;
using std::left;

// Outputs a chunk of information about the simulation, in one line
// time v_sum avg_e_kin avg_e_pot e_tot temp
void output_summary(SimpleMDBox& box, ofstream& ofile) {
    double avg_kinetic = box.kinetic_energy() / box.n_particles;
    double avg_potential = box.potential_energy() / box.n_particles;
    double avg_tot_energy = avg_kinetic + avg_potential;

    ofile << setw(8) << left << box.time
          << setw(13) << left << box.velocity_sum() 
          << setw(12) << left << avg_kinetic
          << setw(12) << left << avg_potential
          << setw(12) << left << avg_tot_energy
          << setw(12) << left << 2 * avg_kinetic / 3
          << "\n";
}

void output_vel_hist(SimpleMDBox& box, ofstream& ofile, double min_val,
                                                        double max_val,
                                                        int num_bins) {
    ofile << box.time << "\t";
    auto histogram = box.compute_velocity_hist(min_val, max_val, num_bins);
    for (auto& bin: histogram) {
        ofile << "(" << bin.first << ", " << bin.second << ")\t";
    }
    ofile << "\n";
}

int main(int argc, char* argv[]) {
    ofstream summ_file("data/summary.dat");
    ofstream vel_hist_file("data/boltzmann.dat");

    const uint N_steps = 10000;
    const vec3 box_dim = {11.0, 11.0, 11.0}; // Want density of ~0.8
    const uint N_particles = 512; // 8x8x8
    const double temp = 1.0;
    const double dt = 0.0005;
    SimpleMDBox box(box_dim, N_particles, temp, dt);
    box.set_langevin(true);

    summ_file << "time    v_sum        avg_e_kin   avg_e_pot   avg_e_tot   temperature\n";
    vel_hist_file << "time\t(bin_left_edge, rel_freq) ...\n";

    for (int n = 0; n < N_steps; ++n) {
        output_summary(box, summ_file);
        if (n % 200 == 0 && n <= 800 || n % 2000 == 0) {
            output_vel_hist(box, vel_hist_file, 0.0, 5.0, 50);
        }
        if (n % 100 == 0) {
            cout << "Step: " << n << "\n";
        }
        box.nvt_step();
    }

    return 0;
}