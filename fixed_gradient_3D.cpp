/*
Fixed chemoattractant concentration. Cells move towords high concentration of varying chemoattracted, the movement is directed,
 but if cells do not sense higher concentration, they move randomly. To change the form of the fixed chemoattractant concentration,
 look at line 82
*/

#include "Aboria.h"
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
#include <iostream>// writing on a text file
#include <fstream>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXf, MatrixXf

//VectorXi func(double diff_conc, double slope, int n_seed) {
int main(){

    //double diff_conc=0.08;
    double slope = 0.0333;
    int n_seed_y = 1;
    int n_seed_z = 5;

    // model parameters

    int length_x = 30;//240;//40;//4; // length of the chemoattractant vector, for fixed domain
    double domain_length = 30; //this variable is for the actual domain length
    double old_length = 30;// this is important for the update of the positions
    const int length_y = 12;//120;//20;//4;
    const int length_z = 12;
    double cell_radius = 0.75;//0.5; // radius of a cell
    const double diameter =
            2 * cell_radius;//2 // diameter in which there have to be no cells, equivalent to size of the cell
    const int N_steps = 800; // number of times the cells move up the gradient
    const size_t N = 7; // initial number of cells
    double l_filo = 4.75;//2; // sensing radius
    double diff_conc = 0.005; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 1; // determines how frequently domain grows (actually not relevant because it will go every timestep)
    int insertion_freq = 5;
    double speed_l = 0.5;//0.05; // speed of a leader cell
    double speed_f = 0.5;//0.08; // speed of a follower cell
    double dettach_prob = 0.5; // probability that a follower cell which is on trail looses the trail


    // distance to the track parameters
    double dist_thres = 1;
    int closest_time;
    int leader_track;


    // domain growth parameters


    /*
    double L_0 = initial_domain_length;//300;
    double a = 0.08/10;
     int L_inf = 50;//1100;//100;//16;//870;
    */

    // for logistic growth with delay

    // correct, but for now use different ones
//    double L_0 = 30; // will have to make this consistent with actual initial length
//    double a = 0.23/60;//0.008;//0.23/10;
//    double L_inf = 86.76;
//    double t_s = 16*60;//4.31*10;
//    double constant = 29.12;

    double L_0 = 30; // will have to make this consistent with actual initial length
    //double a = 0.001;//0.008;//0.23/10;
    double a = 0.23 / 20;
    double L_inf = 86.76;
    //double t_s = 16;//4.31*10;
    double t_s = 15.96 * 20;
    double constant = 29.12;

    double domain_len_der = 0; // initialise derivative of the domain growth function

    // parameters for the dynamics of chemoattractant concentration


    double D = 1; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0; // initialise time, redundant
    double dt = 0.0001; // time step
    double dx = 1; // space step in x direction, double to be consistent with other types
    double dy = 1; // space step in y direction
    double kai = 1 / 100;//0.0001/10; // to 1 /h production rate of chemoattractant


    // parameters for internalisation

    double R = cell_radius;//7.5/10; // \nu m cell radius
    //double lam = 500;//(100)/10; // to 1000 /h chemoattractant internalisation
    // matrix that stores the values of concentration of chemoattractant
    double chemo[length_x][length_y][length_z];



    // generate initial concentration of chemoattractant
    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            for (int k = 0; k < length_z; k++) {
                chemo[i][j][k] =
                        slope * i;//i*i;//log(double(i+1)); // concentration grows linearly/ quadratic/ logistic
            }
        }
    }


    // four columns for x, y, z, u (z is necessaty for paraview)

    // form a matrix which would store x,y,z,u
    MatrixXf chemo_3col = MatrixXf::Zero(length_x * length_y * length_z, 4);
    MatrixXi chemo_3col_ind = MatrixXi::Zero(length_x * length_y * length_z,
                                             3); // need for because that is how paraview accepts data, third dimension is just zeros


    // x, y coord, 1st and 2nd columns respectively
    int k = 0;
    // it has to be 3D for paraview
    while (k < length_x * length_y * length_z) {
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                for (int l = 0; l < length_z; l++) {
                    chemo_3col_ind(k, 0) = i;
                    chemo_3col_ind(k, 1) = j;
                    chemo_3col_ind(k, 2) = l;
                    k += 1;
                }
            }
        }
    }


    // y and x (initially) column
    for (int i = 0; i < length_x * length_y * length_z; i++) {
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);
        chemo_3col(i, 0) = chemo_3col_ind(i, 0);
        chemo_3col(i, 2) = chemo_3col_ind(i, 2);
    }


    // u column
    for (int i = 0; i < length_x * length_y * length_y; i++) {
        chemo_3col(i, 3) = chemo[chemo_3col_ind(i, 0)][chemo_3col_ind(i, 1)][chemo_3col_ind(i, 2)];
    }

    // save data to plot chemoattractant concentration in MATLAB
    ofstream output("matrix_3col.csv");

    output << "x, y, z, u" << "\n" << endl;


    for (int i = 0; i < length_x * length_y * length_z; i++) {
        for (int j = 0; j < 4; j++) {
            output << chemo_3col(i, j) << ", ";
        }
        output << "\n" << endl;
    }

    /*
     * 2D domain with a few randomly placed particles
     */

    /*
     * initial cells of fixed radius
     */

    //ABORIA_VARIABLE(velocity,vdouble2,"velocity")
    ABORIA_VARIABLE(radius, double, "radius")
    typedef Particles<std::tuple<radius>, 3> particle_type;
    //typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
    typedef particle_type::position position;
    particle_type particles(N);
    // uniformly in y
    std::default_random_engine gen;
     gen.seed(n_seed_y);
    std::uniform_real_distribution<double> uniform(2, length_y - 1);

    // uniformly in z
    std::default_random_engine genz;
     genz.seed(n_seed_z);
    std::uniform_real_distribution<double> uniformz(2, length_z - 1);


    /*
         * initialise neighbour search with 2d cuboid domain,
         * periodic in x and y
         */

    // particles.init_neighbour_search(vdouble2(0, 0), vdouble2(5*length_x, length_y), vbool2(false, false));

// only needed here if not compact initialisation
    particles.init_neighbour_search(vdouble3(0, 0, 0), 5 * vdouble3(length_x, length_y, length_z),
                                    vbool3(false, false, false));

    /*
     * random initialisation
     */

    for (int i = 0; i < N; ++i) {
        bool free_position = false;
        get<radius>(particles[i]) = cell_radius;

        get<position>(particles[i]) = vdouble3(cell_radius, uniform(gen),uniformz(genz)); // x=2, uniformly in y


    }



    /*
     * compact initialisation
     */

//    for (int i = 0; i < N; ++i) {
//
//
//        get<radius>(particles[i]) = cell_radius;
//        //get<type>(particles[i]) = 0; // initially all cells are leaders
//
//        //get<position>(p) = vdouble2(cell_radius,(i+1)*diameter); // x=2, uniformly in y
//        get<position>(particles[i]) = vdouble3(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                            0.5 * double(length_y - 1) /
//                                                            double(N), (i + 1) * double(length_y - 1) / double(N) -
//                                                                       0.5 * double(length_y - 1) /
//                                                                       double(N)); // x=2, uniformly in y
//
//
//    }



    // save particles before they move
    vtkWriteGrid("before", 0, particles.get_grid(true));
    vtkWriteGrid("particles", t, particles.get_grid(true));


    // choose a set of random number between 0 and 2*pi
    std::default_random_engine gen1;
    std::uniform_real_distribution<double> uniformpi(0, 2 * M_PI); // angle for theta

    std::default_random_engine genphi;
    std::uniform_real_distribution<double> uniformphi(0, M_PI); // angle for phi



    for (int t = 0; t < N_steps; t++) {
        cout << "domain length " << domain_length << endl;

        // insert new cells at the start of the domain at insertion time (have to think about this insertion time)

        if (t % insertion_freq == 0) {
            bool free_position = false;
            particle_type::value_type p;
            get<radius>(p) = cell_radius;

            get<position>(p) = vdouble3(cell_radius, uniform(gen), uniformz(genz)); // x=2, uniformly in y
            free_position = true;
            /*
             * loop over all neighbouring particles within "dem_diameter" distance
             */
            //particle_type::value_type closest_neighbour;
            for (auto tpl = euclidean_search(particles.get_query(), get<position>(p), diameter); tpl != false; ++tpl) {

                free_position = false;
            }
            //cout << "position from j " << get<id>(closest_neighbour) << endl;
            if (free_position) {
                particles.push_back(p);
            }
        }


        /////////////////////////////////////
        // grow domain


        domain_len_der = 0;

        if (t % freq_growth == 0) {
            //cout << "are you never in " << endl;

            // no time delay
            //domain_length = ((L_inf*exp(a*t))/ (L_inf/L_0 + exp(a*t) - 1) );

            //domain_len_der = ((a*L_inf*exp(a*t))/ (L_inf/L_0 + exp(a*t) - 1) - (a*L_inf*exp(2*a*t))/(L_inf/L_0 + exp(a*t) - 1) );

            // from the paper
            //domain_length = ((L_inf*exp(a*(t-t_s)))/ (L_inf/L_0 + exp(a*(t-t_s)) - 1) ) + constant;

            //domain_len_der = ((a*L_inf*exp(a*(t-t_s)))/ (L_inf/L_0 + exp(a*(t-t_s)) - 1) - (a*L_inf*exp(2*a*(t-t_s)))/(L_inf/L_0 + exp(a*(t-t_s)) - 1) );

            // with time delay and constant to make the initial conditions consistent
            domain_length = ((L_inf * exp(a * (t - t_s))) / (L_inf / L_0 + exp(a * (t - t_s)) - 1)) + constant;

            domain_len_der = ((a * L_inf * exp(a * (t - t_s))) / (L_inf / L_0 + exp(a * (t - t_s)) - 1) -
                              (a * L_inf * exp(2 * a * (t - t_s))) / (L_inf / L_0 + exp(a * (t - t_s)) - 1));


        }

        // update chemoattractant profile

        /*
            Rescale x coordinates properly
        */


        // three columns for x, y, z

        // form a matrix which would store x,y,z, y -same as before


        //cout << "old length " << old_length << endl;
        //cout << "domain length " << domain_length << endl;

        /*
         *  x column, UNIFORM GROWTH
         *
         * */
        for (int i = 0; i < length_x * length_y * length_z; i++) {
            chemo_3col(i, 0) = chemo_3col_ind(i, 0) * (domain_length / length_x);
        }
        cout << "domain length ratio " << domain_length / length_x << endl;



        /*
         *  x column, NON-UNIFORM GROWTH
         *
         * */
        // only the first half grows

        // domain does not grow in the final half
        // the last element is dependent on how much the domain grew
//        for (int j=1;j< length_y+1;j++){
//            chemo_3col(length_x*length_y - j,0) = chemo_3col_ind(length_x*length_y - j,0)*(domain_length / length_x);
//        }
//        // since the second half does not grow, the chemo remains fixed
//        float difference = chemo_3col_ind(length_x*length_y-1,0) - chemo_3col_ind(length_x*length_y-(length_y+1),0);
//        int count = 0.5*length_x;
//        int count_12 = 1; // count, so that the change would be at every twelth position
//        for (int i = 0.5 * length_x * length_y; i < length_x * length_y - length_y; i++) {
//            chemo_3col(i, 0) = chemo_3col(length_x*length_y -1,0) -difference *count;
//            count_12 += 1;
//            cout << " x coord, 2nd half " << chemo_3col(i,0) << endl;
//
//            if (count_12 % length_y == 0){
//                count -= 1;
//            }
//        }
//
//        // the first half grows
//        for (int i = 0; i < 0.5 * length_x * length_y; i++) {
//            chemo_3col(i, 0) = chemo_3col_ind(i, 0) * (chemo_3col(0.5*length_x*length_y,0) / (length_x*0.5));
//            cout << " x coord, 1st half " << chemo_3col(i,0) << endl;
//        }




        // u column

        for (int i = 0; i < length_x * length_y * length_z; i++) {
            chemo_3col(i, 3) = chemo[chemo_3col_ind(i, 0)][chemo_3col_ind(i, 1)][chemo_3col_ind(i, 2)];
        }


        // save data to plot chemoattractant concentration
        ofstream output("matrix_growing_domain" + to_string(t) + ".csv");

        output << "x, y, z, u" << "\n" << endl;


        for (int i = 0; i < length_x * length_y * length_z; i++) {
            for (int j = 0; j < 4; j++) {
                output << chemo_3col(i, j) << ", ";
            }
            output << "\n" << endl;
        }


        // update positions

        /*
         * Uniform growth
         * */

        if (t % freq_growth == 0) {
            for (int i = 0; i < particles.size(); i++) {
                get<position>(particles)[i] *= vdouble3((domain_length / old_length), 1, 1);
            }
            old_length = domain_length;
        }

        /*
         * Non-uniform growth, second part does not grow
         * */

//        if (t % freq_growth == 0) {
//            for (int i = 0; i < particles.size(); i++) {
//                vdouble2 x = get<position>(particles)[i]; // so that I could extract x coord
//                // if in the first part of the domain
//                if (x[0] > 0 && x[0] < old_length - 0.5 * length_x){
//                    get<position>(particles)[i] *= vdouble2((domain_length - length_x*0.5)/(old_length - length_x*0.5),1);//uniform growth in the first part of the domain
//                }
//                //if (x[0] > old_length - 0.5 * length_x && x[0] < old_length){
//                    else{
//                    get<position>(particles)[i] += vdouble2(domain_length - old_length,0);//uniform growth in the first part of the domain
//                }
//
//            }
//            old_length = domain_length;
//        }



        /////////////////////////////////////


        // Update positions based on the gradient




        // SEARCH MULTIPLE TIMES

        // Random ordering

        int no_replacements = particles.size(); // to store the number of cells that have already been picked randomly

        int check_rep = 0; // check for repetitions, 0 no rep, 1 rep

        //for (int j=0; j<N_steps;j++){

        std::default_random_engine gen2;
        gen2.seed(t);
        std::uniform_real_distribution<double> uniform_particles(0, no_replacements); // can only move forward

        VectorXi particle_id = VectorXi::Zero(particles.size());

        for (int i = 0; i < particles.size(); i++) {

            check_rep = 1; // set to 1 to enter the while loop
            while (check_rep == 1) {
                check_rep = 0; // it is initially zero and then will be changed to 1 if it is equivalent to others
                particle_id(i) = uniform_particles(gen2);


                for (int j = 0; j < i; j++) {
                    if (particle_id(i) == particle_id(j)) { check_rep = 1; }
                }
                cout << "particle id before " << particle_id(i) << endl;
            }
            //cout << "ids " << particle_id(i) << endl;
        }




        // ordered list

//        for (int i = 0; i < particles.size(); i++) {
//
//                particle_id(i) = i;
//
//        }


//            // find a particle with particle_id
//
//            for (int i = 0; i < particles.size(); i++) {
//                if (get<id>(particles[i]) == particle_id){
//                    int current_id = particle_id;
//                }
//            }



        // pick a cell randomly
        //int particle_id(j) = 0;

        for (int j = 0; j < particles.size(); j++) {
            cout << " consider this id " << particle_id(j) << endl;


//        }
//
//
//        cout << "the number of particles in the system " << particles.size() << endl;
//        for (int j = 0; j < particles.size(); j++ ) {
//            cout << " consider this id " << particle_id(j) << endl;
//            particle_id(j) = particle_id(j);
//            for (int i_bef = 0; i_bef < particles.size(); i_bef++) {
//                if (particle_id(j) == get<id>(particles[i_bef])){
//                    particle_id(j) = i_bef;
//                }
//            }


            //for (int i = 0; i < particles.size(); i++) {

//            cout << "print which particle id I am looping over now " << get<id>(particles[i]) << endl;
//            cout << "print particle size " << particles.size() << endl;
//

            vdouble3 x;
            x = get<position>(particles[particle_id(j)]);
            //cout << "particles.size " << i << endl;
            //cout << "print id " << get<id>(particles[i]) << endl;
            //cout << "x coord " << x[0] << endl;



            /*
             * x_in variable will correspond to the coordinates on the non-updated domain (same as initial)
             * */

            // Uniform domain growth

            double x_in; // x coordinate in initial domain length

            x_in = (length_x / domain_length) * x[0];

            // Non-uniform domain growth



            // Non-uniform domain growth, onl first half grows
            // if in the first part of the domain
//            if (x[0] > 0 && x[0] < domain_length - 0.5 * length_x){
//                x_in = x[0]*((length_x - length_x*0.5)/(domain_length - length_x*0.5));//uniform growth in the first part of the domain
//            }
//            if (x[0] > domain_length - 0.5 * length_x && x[0] < domain_length){
//                x_in = x[0]- (domain_length - length_x);//uniform growth in the first part of the domain
//            }





            // create an array to store random directions
            std::array<double, 3> random_angle;
            array<double, 3> random_angle_phi;
//            std::array<int, 3> sign_x;
//            std::array<int, 3> sign_y;
            for (int i = 0; i < 3; i++) {

                double random_angle_tem = uniformpi(gen1);
                double random_angle_tem_phi = uniformphi(genphi);
                // int sign_x_tem, sign_y_tem;

                while (round((x_in + cos(random_angle_tem) * sin(random_angle_tem_phi) * l_filo)) < 0 ||
                       round((x_in + cos(random_angle_tem) * sin(random_angle_tem_phi) * l_filo)) >
                       length_x - 1 || round(x[1] + sin(random_angle_tem) * sin(random_angle_tem_phi) * l_filo) < 0 ||
                       round(x[1] + cos(random_angle_tem) * sin(random_angle_tem_phi) * l_filo) > length_y - 1 ||
                       round(x[2] + cos(random_angle_tem_phi) * l_filo) < 0 ||
                       round(x[2] + cos(random_angle_tem_phi) * l_filo) >
                       length_z - 1) {
                    random_angle_tem = uniformpi(gen1);
                    random_angle_tem_phi = uniformphi(genphi);

//                    if (sin(random_angle_tem) < 0) {
//                        sign_x_tem = -1;
//                    } else { sign_x_tem = 1; }
//
//                    if (cos(random_angle_tem) < 0) {
//                        sign_y_tem = -1;
//                    } else { sign_y_tem = 1; }

                }

                random_angle[i] = random_angle_tem;
                random_angle_phi[i] = random_angle_tem_phi;
                cout << "theta " << random_angle[i] << endl;
                cout << "phi " << random_angle_phi[i] << endl;

//                sign_x[i] = sign_x_tem;
//                sign_y[i] = sign_y_tem;


            }


            // choose which direction to move


            cout << "chemo coef x before " << int(round(x_in)) << endl;
            cout << "chemo coef x after " << int(round(x_in + cos(random_angle[0]) * sin(random_angle_phi[0]) * l_filo)) << endl;

            cout << "chemo coef y before " << int(round(x[1])) << endl;
            cout << "chemo coef y after " << int(round(x[1]+ cos(random_angle[0]) * sin(random_angle_phi[0]) * l_filo)) << endl;

            cout << "chemo coef z before " << int(round(x[2])) << endl;
            cout << "chemo coef z after " << int(round(x[2] + cos(random_angle[0]) * sin(random_angle_phi[0]) * l_filo)) << endl;




            // store variables for concentration at new locations
            double old_chemo = chemo[int(round(x_in))][int(round(x)[1])][int(round(x)[2])];
            //double new_chemo_1 = chemo(round((x[0] * (length_x / domain_length) + sin(random_angle[0]) + sign_x[0] * l_filo)),
            //round(x[1] + cos(random_angle[0]) + sign_y[0] * l_filo));
            double new_chemo_1 = chemo[int(round(x_in + cos(random_angle[0]) * sin(random_angle_phi[0]) * l_filo))][int(
                    round(
                            x[1] + sin(random_angle[0]) * sin(random_angle_phi[0]) * l_filo))][int(round(
                    x[2] + cos(random_angle_phi[0]) * l_filo))];

            //double new_chemo_2 = chemo(round((x[0] * (length_x / domain_length) + sin(random_angle[1]) + sign_x[1] * l_filo)),
            //round(x[1] + cos(random_angle[1]) + sign_y[1] * l_filo));
            double new_chemo_2 = chemo[int(
                    round(x_in + cos(random_angle[1]) * sin(random_angle_phi[1]) * l_filo))][int(
                    round(x[1] + sin(random_angle[1]) * sin(random_angle_phi[1]) * l_filo))][int(
                    round(x[2] + cos(random_angle_phi[1]) * l_filo))];


            cout << "new_chemo_2 " << new_chemo_2 << endl;
            cout << "new_chemo_1 " << new_chemo_1 << endl;
            cout << "old " << old_chemo << endl;



            //cout << "stops here " << endl;
            // if first direction greater, second smaller
            //absolute
            //else if (new_chemo_1 - old_chemo > diff_conc && new_chemo_2 - old_chemo < diff_conc){

            //relative

            if ((new_chemo_1 - old_chemo) / sqrt(old_chemo) > diff_conc &&
                (new_chemo_2 - old_chemo) / sqrt(old_chemo) < diff_conc) {


                x += speed_l * vdouble3(cos(random_angle[0]) * sin(random_angle_phi[0]),
                                        sin(random_angle[0]) * sin(random_angle_phi[0]), cos(random_angle_phi[0]));

                //cout << "print id " << id_[x] << endl;
                x_in = (length_x / domain_length) * x[0];


                //cout << "Position "<< x << endl;
                int count_position = 0;
                bool free_position = true; // check if the neighbouring position is free

                // if this loop is entered, it means that there is another cell where I want to move
                for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; k++) {


                    if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                    //break;
                }


                //cout << "print position " << count_position << endl;
                // check that the position they want to move to is free and not out of bounds
                // check that the position they want to move to is free and not out of bounds
                if (free_position && round(x_in) >= 0 &&
                    round(x_in) <= length_x - 1 && round(x[1]) >= 0 &&
                    round(x[1]) <= length_y - 1 && round(x[2]) >= 0 &&
                    round(x[2]) <= length_y - 1) {
                    get<position>(particles)[particle_id(j)] +=
                            speed_l * vdouble3(cos(random_angle[0]) * sin(random_angle_phi[0]),
                                               sin(random_angle[0]) * sin(random_angle_phi[0]),
                                               cos(random_angle_phi[0])); // update if nothing is in the next position
                }

            }
                // if first smaller, second bigger

                //absolute
                // else if (new_chemo_1 - old_chemo < diff_conc && new_chemo_2 - old_chemo > diff_conc){

                //relative
            else if ((new_chemo_1 - old_chemo) / sqrt(old_chemo) < diff_conc &&
                     (new_chemo_2 - old_chemo) / sqrt(old_chemo) > diff_conc) {

                x += speed_l * vdouble3(cos(random_angle[1]) * sin(random_angle_phi[1]),
                                        sin(random_angle[1]) * sin(random_angle_phi[1]), cos(random_angle_phi[1]));
                //cout << "print id " << id_[x] << endl;
                x_in = (length_x / domain_length) * x[0];

                //cout << "Position "<< x << endl;
                int count_position = 0;
                bool free_position = true; // check if the neighbouring position is free

                // if this loop is entered, it means that there is another cell where I want to move
                for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; k++) {


                    if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                    //break;
                }


                //cout << "print position " << count_position << endl;
                // check that the position they want to move to is free and not out of bounds
                // check that the position they want to move to is free and not out of bounds
                if (free_position && round(x_in) >= 0 &&
                    round(x_in) <= length_x - 1 && round(x[1]) >= 0 &&
                    round(x[1]) <= length_y - 1 && round(x[2]) >= 0 &&
                    round(x[2]) <= length_y - 1) {
                    get<position>(particles)[particle_id(j)] +=
                            speed_l * vdouble3(cos(random_angle[1]) * sin(random_angle_phi[1]),
                                               sin(random_angle[1]) * sin(random_angle_phi[1]),
                                               cos(random_angle_phi[1])); // update if nothing is in the next position
                }
                //break;
            }
                // if both greater choose the bigger one

                // absolute
                //else if (new_chemo_1 - old_chemo > diff_conc && new_chemo_2 - old_chemo > diff_conc){

                //relative
            else if ((new_chemo_1 - old_chemo) / sqrt(old_chemo) > diff_conc &&
                     (new_chemo_2 - old_chemo) / sqrt(old_chemo) > diff_conc) {


                // if first is greater than the second
                if (new_chemo_1 > new_chemo_2) {

                    x += speed_l * vdouble3(cos(random_angle[0]) * sin(random_angle_phi[0]),
                                            sin(random_angle[0]) * sin(random_angle_phi[0]), cos(random_angle_phi[0]));

                    //cout << "print id " << id_[x] << endl;
                    x_in = (length_x / domain_length) * x[0];

                    //cout << "Position "<< x << endl;
                    int count_position = 0;
                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move
                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; k++) {


                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                        }
                        //break;
                    }
                    //cout << "print position " << count_position << endl;
                    // check that the position they want to move to is free and not out of bounds
                    // check that the position they want to move to is free and not out of bounds
                    if (free_position && round(x_in) >= 0 &&
                        round(x_in) <= length_x - 1 && round(x[1]) >= 0 &&
                        round(x[1]) <= length_y - 1 && round(x[2]) >= 0 &&
                        round(x[2]) <= length_y - 1) {
                        get<position>(particles)[particle_id(j)] +=
                                speed_l * vdouble3(cos(random_angle[0]) * sin(random_angle_phi[0]),
                                                   sin(random_angle[0]) * sin(random_angle_phi[0]),
                                                   cos(random_angle_phi[0])); // update if nothing is in the next position
                    }

                }
                    // if second is greater than the first
                else if (new_chemo_1 < new_chemo_2) {

                    x += speed_l * vdouble3(cos(random_angle[1]) * sin(random_angle_phi[1]),
                                            sin(random_angle[1]) * sin(random_angle_phi[1]), cos(random_angle_phi[1]));
                    //cout << "print id " << id_[x] << endl;
                    x_in = (length_x / domain_length) * x[0];

                    //cout << "Position "<< x << endl;
                    int count_position = 0;
                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move
                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; k++) {


                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                        }
                        //break;
                    }

                    //cout << "print position " << count_position << endl;
                    // check that the position they want to move to is free and not out of bounds
                    // check that the position they want to move to is free and not out of bounds
                    if (free_position && round(x_in) >= 0 &&
                        round(x_in) <= length_x - 1 && round(x[1]) >= 0 &&
                        round(x[1]) <= length_y - 1 && round(x[2]) >= 0 &&
                        round(x[2]) <= length_y - 1) {
                        get<position>(particles)[particle_id(j)] +=
                                speed_l * vdouble3(cos(random_angle[1]) * sin(random_angle_phi[1]),
                                                   sin(random_angle[1]) * sin(random_angle_phi[1]),
                                                   cos(random_angle_phi[1])); // update if nothing is in the next position
                    }

                }

            }



            //if both smaller, move random direction
            //absolute
            //if (new_chemo_1 - old_chemo < diff_conc && new_chemo_2 - old_chemo < diff_conc) {


            // relative
            else{

                cout << "is it ever here " << endl;
                x += speed_l * vdouble3(cos(random_angle[2]) * sin(random_angle_phi[2]),
                                        sin(random_angle[2]) * sin(random_angle_phi[2]), cos(random_angle_phi[2]));
                //cout << "print id " << id_[x] << endl;
                x_in = (length_x / domain_length) * x[0];


                //cout << "Position "<< x << endl;
                int count_position = 0;
                bool free_position = true; // check if the neighbouring position is free

                // if this loop is entered, it means that there is another cell where I want to move
                for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; k++) {


                    if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                    //break;
                }


                //cout << "print position " << count_position << endl;

                cout << " position before the update " << get<position>(particles)[particle_id(j)] << endl;
                // check that the position they want to move to is free and not out of bounds
                if (free_position && round(x_in) >= 0 &&
                    round(x_in) <= length_x - 1 && round(x[1]) >= 0 &&
                    round(x[1]) <= length_y - 1 && round(x[2]) >= 0 &&
                    round(x[2]) <= length_y - 1) {
                    cout << " and here " << endl;
                    get<position>(particles)[particle_id(j)] +=
                            speed_l * vdouble3(cos(random_angle[2]) * sin(random_angle_phi[2]),
                                               sin(random_angle[2]) * sin(random_angle_phi[2]),
                                               cos(random_angle_phi[2])); // update if nothing is in the next position

                    cout << " position after the update " << get<position>(particles)[particle_id(j)] << endl;

                }

            }


            //cout << "sin of the angle " << sin(random_angle(rand_num_count)) << endl;
            //cout << "cos of the angle " << cos(random_angle(rand_num_count)) << endl;


            // check if the x coordinates are not out of the domain, if they are, ignore that step

            //cout << "chemo coord x before the test " << round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_radius) << endl;
            //cout << "chemo coord y before the test " << round(x[0]+cos(random_angle(rand_num_count))+sign_y*cell_radius) << endl;

            //rand_num_count += 3; // update random number count

            //}// go through all the particles

        } // go through all ids in a random vector, thus all the particles



        particles.update_positions();

        //particles.update_positions();
        // save particles after they move
        /*
     * on every i/o step write particle container to a vtk
     * unstructured grid file
     */
        //cout << "." << flush;
#ifdef HAVE_VTK
        vtkWriteGrid("particles", t + 1, particles.get_grid(true));
#endif

        cout << " end of a time step " << endl;
    }//all time steps




    // count the number of particles that are at the last 10% of the domain

    /*
    double last_10_domain = domain_length - domain_length*domain_fraction;
    cout << "last 10 percent " << last_10_domain << endl;
    int number_of_cells = 0; // number of cells in the last 10% of the domain

    for (int p =0; p< particles.size(); p++){
        vdouble2 x = get<position>(particles[p]);
        if(x[0] > last_10_domain )
        number_of_cells += 1;
    }


    double proportion_cells_last = double(number_of_cells)/double(particles.size());

    cout << "proportion of cells " << proportion_cells_last << endl;

    return proportion_cells_last;
    */


    /*
    return the furthest distance travelled by the cells
     */

//    double furthest_distance = 0 ;
//    vdouble2 dist; // variable for positions
//
//    for (int i = 0; i < particles.size(); i++) {
//
//        dist = get<position>(particles[i]);
//
//        if (furthest_distance < dist[0]){
//            furthest_distance = dist[0];
//        }
//
//    }
//
//    return furthest_distance;

    /*
     * return the density of cells in each of the fifth of the domain
     */
//    const int domain_partition = int(domain_length / double(5));; // number of intervalas of 50 \mu m
//
//    VectorXi proportions = VectorXi::Zero(domain_partition); // integer with number of cells in particular part
//    //array<double, domain_partition> proportions;
//
//
//    double one_part = domain_length / double(domain_partition);
//
//    cout << "one part of the domain " << one_part << endl;
//
//    for (int i = 0; i < domain_partition; i++) {
//
//        for (int j = 0; j < particles.size(); j++) {
//            vdouble3 x = get<position>(particles[j]);
//            //cout<< "domain partition " << i*one_part << endl;
//            //cout << "x coordinate " << x[0] << endl;
//            if (i * one_part < x[0] && x[0] < (i + 1) * one_part) {
//                proportions(i) += 1;
//            }
//        }
//
//    }
//
//
//    // for loops to count the number of cells in each of the fifth of the domain
//
//    return proportions;


}


// parameter analysis

///*
// * main for futhest distance
// */
//
//int main(){
//
//    const int number_parameters = 1; // parameter range
//    const int sim_num = 1;
//
//    MatrixXf all_distances = MatrixXf::Zero(number_parameters,sim_num); //matrix over which I am going to average
//
////n would correspond to different seeds
//    for (int n = 0; n < sim_num; n++) {
//
//
//        // define parameters that I will change
//        //VectorXf slope, threshold;
//        array<double, number_parameters> threshold;
//        array<double, 1> slope;
//        //array<double,number_parameters,number_parameters>;
//
//        MatrixXf furthest_distance = MatrixXf::Zero(number_parameters,1);
//        //VectorXf furthest_distance = VectorXf::Zero(number_parameters);
//
//
//        for (int i = 0; i < number_parameters; i++) {
//            threshold[i] = 0.05; // just an example
//            //threshold[i] = 0.005 * (i + 1);// 0.01;
//            //cout << "slope " << slope[i] << endl;
//
//        }
//
//        slope[0] = 0.0333; //linear growth ax
//        //slope[0] = 0.0011; // quadratic growth ax^2
//        //slope[0] = 0.2912; // logistic growth
//
//
//        // VectorXi numbers = func(0.005, slope[0], 2);
//
//
//        for (int i = 0; i < number_parameters; i++) {
//
//            //for (int j = 0; j < 1; j++) {
//
//            furthest_distance(i, 0) = func(threshold[i], slope[0], n);
//            cout << "number of parameters investigated " << i << endl;
//
//            //}
//            all_distances(i,n) = furthest_distance(i,0);
//        }
//
//
//        // save data to plot chemoattractant concentration
////        ofstream output("furthest_distance_matrix.csv");
////
////        MatrixXf furthest_distance_3col(number_parameters * number_parameters, 4), furthest_distance_3col_ind(number_parameters * number_parameters,
////                                                                    2); // need for because that is how paraview accepts data, third dimension is just zeros
////
////
////
////        // x, y coord, 1st and 2nd columns respectively
////        int k = 0;
////        // it has to be 3D for paraview
////        while (k < number_parameters * number_parameters) {
////            for (int i = 0; i < number_parameters; i++) {
////                for (int j = 0; j < number_parameters; j++) {
////                    furthest_distance_3col_ind(k, 0) = i;
////                    furthest_distance_3col_ind(k, 1) = j;
////                    furthest_distance_3col(k, 2) = 0;
////                    k += 1;
////                }
////            }
////        }
////
////
////        // y and x (initially) column
////        for (int i = 0; i < number_parameters * number_parameters; i++) {
////            furthest_distance_3col(i, 1) = furthest_distance_3col_ind(i, 1);
////            furthest_distance_3col(i, 0) = furthest_distance_3col_ind(i, 0);
////        }
////
////
////        // u column
////        for (int i = 0; i < number_parameters * number_parameters; i++) {
////            furthest_distance_3col(i, 3) = furthest_distance(furthest_distance_3col_ind(i, 0), furthest_distance_3col_ind(i, 1));
////        }
////
////        output << "x, y, z, u" << "\n" << endl;
////
////
////        for (int i = 0; i < number_parameters * number_parameters; i++) {
////            for (int j = 0; j < 4; j++) {
////                output << furthest_distance_3col(i, j) << ", ";
////            }
////            output << "\n" << endl;
////        }
////
//
//
//
//
//        // This is what I am using for MATLAB
//        ofstream output2("furthest_distance_matrix_matlab.csv");
//
//        for (int i = 0; i < number_parameters; i++) {
//
//            for (int j = 0; j < 1; j++) {
//
//                output2 << furthest_distance(i, j) << ", ";
//
//            }
//            output2 << "\n" << endl;
//        }
//
//    }
//
//
//
//    ofstream output3("simulations_simple.csv");
//
//    for (int i = 0; i < number_parameters; i++) {
//
//        for (int j = 0; j < sim_num; j++) {
//
//            output3 << all_distances(i, j) << ", ";
//
//        }
//        output3 << "\n" << endl;
//    }
//
//
//}

/*
 * main for proportions in different sections
 */


// parameter analysis
//int main() {
//
//    const int number_parameters = 1; // parameter range
//    const int sim_num = 1;
//
//    VectorXi vector_check_length = func(0.005, 0.1, 2); //just to know what the length is
//
//    int num_parts = vector_check_length.size(); // number of parts that I partition my domain
//
//    MatrixXf sum_of_all = MatrixXf::Zero(num_parts, number_parameters); // sum of the values over all simulations
//
////n would correspond to different seeds
//    for (int n = 0; n < sim_num; n++) {
//
//
//        // define parameters that I will change
//        //VectorXf slope, threshold;
//        array<double, number_parameters> threshold;
//        array<double, 1> slope;
//        //array<double,number_parameters,number_parameters>;
//
//        MatrixXf furthest_distance = MatrixXf::Zero(number_parameters, 1);
//        //VectorXf furthest_distance = VectorXf::Zero(number_parameters);
//
//
//        for (int i = 0; i < number_parameters; i++) {
//            threshold[i] = 0.5;
//            //threshold[i] = 0.005 * (i + 1);// 0.01;
//            //cout << "slope " << slope[i] << endl;
//
//        }
//
//        slope[0] = 0.0333;//0.0175; //linear growth ax
//        //slope[0] = 0.0011; // quadratic growth ax^2
//        //slope[0] = 0.2912; // logistic growth
//
//
//        MatrixXi numbers = MatrixXi::Zero(num_parts,
//                                          number_parameters); // can't initialise because do not know the size
//
//
//        for (int i = 0; i < number_parameters; i++) {
//
//            //for (int j = 0; j < 1; j++) {
//
//            numbers.block(0, i, num_parts, 1) = func(threshold[i], slope[0], n);
//
//            //}
//        }
//
//
//        // This is what I am using for MATLAB
//        ofstream output2("numbers_matrix_matlab.csv");
//
//        for (int i = 0; i < numbers.rows(); i++) {
//
//            for (int j = 0; j < numbers.cols(); j++) {
//
//                output2 << numbers(i, j) << ", ";
//
//                sum_of_all(i, j) += numbers(i, j);
//
//            }
//            output2 << "\n" << endl;
//        }
//
//    }
//
//    /*
//    * will store everything in one matrix, the entries will be summed over all simulations
//    */
//
//    ofstream output3("simulations_domain_partition_simple.csv");
//
//    for (int i = 0; i < num_parts; i++) {
//
//        for (int j = 0; j < number_parameters; j++) {
//
//            output3 << sum_of_all(i, j) << ", ";
//
//        }
//        output3 << "\n" << endl;
//    }
//
//
//}