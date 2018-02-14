#include "../include/simpleMD.hpp"
#include "../include/tclap/CmdLine.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <chrono>

using std::cout;
using std::cerr;
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

void output_hist(Histogram& hist, ofstream& ofile) {
    for (const auto& bin: hist) {
        ofile << "(" << bin.left_edge << ", " << bin.count << ")\t";
    }
    ofile << "\n";
}

int main(int argc, char* argv[]) {
    
    struct cmdargs {
        bool langevin;
        unsigned int N_steps;
        double box_side;
        unsigned int N_particles;
        double temp;
        double dt;
        std::string ghist_filename;
        std::string summary_filename;
        std::string vhist_filename;
    } args;

    try {
        /* All of this is command line argument parsing, I'm using 
         *   the TCLAP library, which you can find here: 
         *   http://tclap.sourceforge.net/
         * There's no need to install it, though, as I included the 
         *   necessary files in the
         *   project.
         */

        TCLAP::CmdLine cmd("Simple 3D MD Simulation", ' ', "0.1");
        TCLAP::SwitchArg langevin_arg("l", "langevin", 
                    "Enables langevin dynamics.", cmd, false);
        TCLAP::ValueArg<unsigned int>  N_steps_arg("n", "nsteps", 
                    "The number of steps to perform.", false, 5000, 
                    "unsigned int", cmd);
        TCLAP::ValueArg<double> box_side_arg("b", "boxside", 
                    "The length of the sides of the box.", false, 9.0, 
                    "double", cmd);
        TCLAP::ValueArg<unsigned int> N_particles_arg("N", "nparticles", 
                    "The number of particles in the system.", false, 512, 
                    "unsigned int", cmd);
        TCLAP::ValueArg<double> temp_arg("T", "temp", 
                    "The temperature of the system.", false, 1.0, 
                    "double", cmd);
        TCLAP::ValueArg<double> dt_arg("t", "step", "The step size in time.", 
                    false, 0.0005, "double", cmd);
        TCLAP::ValueArg<std::string> ghist_arg("g", "ghist", 
                    "The name of the file to save the g histogram to, in data/", 
                    false, "ghist.dat", "string", cmd);
        TCLAP::ValueArg<std::string> sum_arg("s", "summary",
                    "The name of the file to save the summary to, in data/",
                    false, "summary.dat", "string", cmd);
        TCLAP::ValueArg<std::string> vhist_arg("v", "vhist",
                    "The name of the file to save the vel histogram to, in data/",
                    false, "boltzmann.dat", "string", cmd);

        cmd.parse(argc, argv);

        args.langevin = langevin_arg.getValue();
        args.N_steps = N_steps_arg.getValue();
        args.box_side = box_side_arg.getValue();
        args.N_particles = N_particles_arg.getValue();
        args.temp = temp_arg.getValue();
        args.dt = dt_arg.getValue();
        args.ghist_filename = ghist_arg.getValue();
        args.summary_filename = sum_arg.getValue();
        args.vhist_filename = vhist_arg.getValue();
    } catch (TCLAP::ArgException& e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
        return 1;
    }

    const vec3 box_dim = {args.box_side, args.box_side, args.box_side};
    SimpleMDBox box(box_dim, args.N_particles, args.temp, args.dt);
    box.set_langevin(args.langevin);

    ofstream summ_file("data/" + args.summary_filename);
    ofstream vel_hist_file("data/" + args.vhist_filename);
    ofstream g_hist_file("data/" + args.ghist_filename);

    summ_file << "time    v_sum        avg_e_kin   avg_e_pot   avg_e_tot   temperature\n";
    vel_hist_file << "time\t(bin_left_edge, rel_freq) ...\n";
    g_hist_file << "time\t(bin_left_edge, rel_freq) ...\n";

    auto t1 = std::chrono::high_resolution_clock::now();
    for (int n = 0; n < args.N_steps; ++n) {
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
    auto t2 = std::chrono::high_resolution_clock::now();
    cout << "Finished in " 
         << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() 
             / 1000.0
         << " seconds\n";

    // We have to scale the g hist specially
    const double bw = box.hists["g"].bin_width;
    const double V = std::pow(args.box_side, 3);
    const int N = box.n_particles;
    constexpr double PI = std::atan(1) * 4;
    auto g_scale = [bw, V, N, PI] (double r) {
        return V / (2.0 * PI * bw * std::pow(N * r, 2));
    };

    Histogram vhist = box.hists["v"];
    vhist.normalize();
    Histogram ghist = box.hists["g"];
    ghist.scale(g_scale);
    ghist.normalize();

    output_hist(vhist, vel_hist_file);
    output_hist(ghist, g_hist_file);

    return 0;
}