#include "simpleMD.hpp"
#include <iostream>

using std::cout;
using std::ofstream;

// Outputs a chunk of information about the simulation, in one line
// time v_sum avg_e_kin avg_e_pot e_tot 
void output_summary(SimpleMDBox& box, ofstream& ofile) {
    double avg_kinetic = box.kinetic_energy() / box.n_particles;
    double avg_potential = box.potential_energy() / box.n_particles;
    double total_energy = box.n_particles * (avg_kinetic + avg_potential);

    ofile << box.time << "\t" << box.velocity_sum() 
                      << "\t" << avg_kinetic
                      << "\t" << avg_potential
                      << "\t" << total_energy
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
    ofstream summ_file("summary.dat");
    ofstream vel_hist_file("boltzmann.dat");

    const uint N_steps = 1000;
    const vec3 box_dim = {11.0, 11.0, 11.0}; // Want density of ~0.8
    const uint N_particles = 1000;
    const double temp = 1.0;
    const double dt = 0.0001;
    SimpleMDBox box(box_dim, N_particles, temp, dt);

    summ_file << "time\tv_sum\tavg_e_kin\tavg_e_pot\te_tot\n";
    vel_hist_file << "time\t(bin_left_edge, rel_freq) ...\n";

    for (int n = 0; n < N_steps; ++n) {
        output_summary(box, summ_file);
        if (n % 200 == 0 && n <= 1000) {
            output_vel_hist(box, vel_hist_file, 0.0, 2.0, 50);
        }
        if (n % 100 == 0) {
            cout << "Step: " << n << "\n";
        }
        box.nvt_step();
    }

    return 0;
}