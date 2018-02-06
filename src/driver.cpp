#include "../include/simpleMD.hpp"
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

void output_hist(SimpleMDBox& box, std::string name, ofstream& ofile) {
    ofile << box.time << "\t";
    auto histogram = box.hists.get_normalized_hist(name);
    for (auto& bin: histogram) {
        ofile << "(" << bin.left_edge << ", " << bin.count << ")\t";
    }
    ofile << "\n";
}

int main(int argc, char* argv[]) {
    ofstream summ_file("data/summary.dat");
    ofstream vel_hist_file("data/boltzmann.dat");
    ofstream g_hist_file("data/ghist.dat");

    const unsigned int N_steps = 10000;
    const vec3 box_dim = {9.0, 9.0, 9.0}; // Want density of ~0.8
    const unsigned int N_particles = 512; // 8x8x8
    const double temp = 1.0;
    const double dt = 0.0005;
    SimpleMDBox box(box_dim, N_particles, temp, dt);
    box.set_langevin(true);

    summ_file << "time    v_sum        avg_e_kin   avg_e_pot   avg_e_tot   temperature\n";
    vel_hist_file << "time\t(bin_left_edge, rel_freq) ...\n";
    g_hist_file << "time\t(bin_left_edge, rel_freq) ...\n";

    for (int n = 0; n < N_steps; ++n) {
        output_summary(box, summ_file);
        if (n > 3000 && n % 100 == 0) {
            box.hists.update_hist("v");
            box.hists.update_hist("g");
        }
        if (n % 100 == 0) {
            cout << "Step: " << n << "\n";
        }
        box.nvt_step();
    }

    output_hist(box, "v", vel_hist_file);
    output_hist(box, "g", g_hist_file);

    return 0;
}