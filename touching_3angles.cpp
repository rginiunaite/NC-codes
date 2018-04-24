
/*
cells move towords high concentration of varying chemoattracted, the movement is directed,
 cells only move if they sense higher concentration. Leaders and followers separately, cells move in chains, physical
 forcing making them move all together if they are close enough with a leader in front
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

//VectorXi proportions(double diff_conc, int n_seed) {
int main(){


    // model parameters

    int length_x = 30;//240;//40;//4; // length of the chemoattractant vector, for fixed domain
    double domain_length = 30; //this variable is for the actual domain length
    double old_length = 30;// this is important for the update of the positions
    const int length_y = 12;//120;//20;//4;
    double cell_radius = 0.75;//0.5; // radius of a cell
    const double diameter =
            2 * cell_radius;//2 // diameter in which there have to be no cells, equivalent to size of the cell
    const int N_steps = 500; // number of times the cells move up the gradient
    const size_t N = 5; // initial number of cells
    double l_filo_y = 2.75;//2; // sensing radius
    double l_filo_x = 2.75; // will have to rescale filopodia when domain grows
    double l_filo_max = l_filo_y * 2;
    double diff_conc = 0.05; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 1; // determines how frequently domain grows (actually not relevant because it will go every timestep)
    int insertion_freq = 1;
    double speed_l = 0.2;//0.05; // speed of a leader cell
    double speed_f = 0.2;//0.08; // speed of a follower cell
    double dettach_prob = 0.5; // probability that a follower cell which is on trail looses the trail
    double chemo_leader = 0.9; //0.5; // phenotypic switching happens when the concentration of chemoattractant is higher than this (presentation video 0.95), no phenotypic switching
    double eps = 1; // for phenotypic switching, the distance has to be that much higher
    const int filo_number = 2;


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

//    double L_0 = 30; // will have to make this consistent with actual initial length
//    double a = 0.001;//0.008;//0.23/10;
//    double L_inf = 86.76;
//    double t_s = 16;//4.31*10;
//    double constant = 29.12;

    double L_0 = 30;
    double a = 0.23 / 60;
    double L_inf = 86.76;
    double t_s = 16 * 60;
    double constant = 29.12;

    double domain_len_der = 0; // initialise derivative of the domain growth function

    // parameters for the dynamics of chemoattractant concentration


    double D = 1; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0; // initialise time, redundant
    double dt = 0.00001; // time step
    double dx = 1; // space step in x direction, double to be consistent with other types
    double dy = 1; // space step in y direction
    double kai = 1 / 100;//0.0001/10; // to 1 /h production rate of chemoattractant


    // parameters for internalisation

    double R = cell_radius;//7.5/10; // \nu m cell radius
    double lam = 50;//(100)/10; // to 1000 /h chemoattractant internalisation


    // matrix that stores the values of concentration of chemoattractant
    MatrixXf chemo = MatrixXf::Zero(length_x, length_y);
    MatrixXf chemo_new = MatrixXf::Zero(length_x, length_y);

    // initialise internalisation matrix
    MatrixXf intern = MatrixXf::Zero(length_x, length_y);

    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = 1; // uniform concentration initially
            chemo_new(i, j) = 1; // this is for later updates
        }    // generate initial chemoattractant concentration

    }


    // four columns for x, y, z, u (z is necessary for paraview)

    // form a matrix which would store x,y,z,u

    MatrixXf chemo_3col(length_x * length_y, 4), chemo_3col_ind(length_x * length_y,
                                                                2); // need for because that is how paraview accepts data, third dimension is just zeros

    // x, y coord, 1st and 2nd columns respectively
    int k = 0;
    // it has to be 3D for paraview
    while (k < length_x * length_y) {
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                chemo_3col_ind(k, 0) = i;
                chemo_3col_ind(k, 1) = j;
                chemo_3col(k, 2) = 0;
                k += 1;
            }
        }
    }


    // y and x (initially) column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);
        chemo_3col(i, 0) = chemo_3col_ind(i, 0);
    }



    /*
     * 2D domain with a few randomly placed particles
     */

    /*
     * initial cells of fixed radius
     */

    //ABORIA_VARIABLE(velocity,vdouble2,"velocity")
    ABORIA_VARIABLE(radius, double, "radius")
    ABORIA_VARIABLE(direction, vdouble2, "direction")// stores the direction a particle moved
    ABORIA_VARIABLE(attached_to_id, int, "attached_to_id")
    ABORIA_VARIABLE(type, int, "type") // 0 if a cell is a leader, 1 if follower
    ABORIA_VARIABLE(chain_type, int, "chain_type") // 0 if a cell is a leader, 1 if follower
    ABORIA_VARIABLE(chain, double, "chain") // stores whether attached to a leader or follower
    // stores the distance to the closest neighbour, if less than thresold
    typedef Particles<std::tuple<radius, type, attached_to_id, direction, chain, chain_type>, 2> particle_type; // 2 stands for dimension


    /*
     * if attached to a leader attached_at_time_step = 2;
     * if attached to a follower attached_at_time_step = 1;
     * if dettached attached_at_time_step = 0;
     * */
    //typedef Particles<std::tuple<radius, attached_to, >,2> attached_follower_type;
    //typedef Particles<std::tuple<radius,type>,2> detached_follower_type;

    //typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
    typedef particle_type::position position;

    particle_type particles(N);



    /*for (int i=0; i<N; ++i) {
        array<vdouble2, N_steps> track_distance[i];
    }*/



    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform(2, length_y - 1);

    /*
     * initialise neighbour search with 2d cuboid domain,
     * periodic in x and y
     */

    particles.init_neighbour_search(vdouble2(0, 0), 5 * vdouble2(length_x, length_y), vbool2(false, false));


    /*
     * compact initialisation
     */

    for (int i = 0; i < N; ++i) {


        get<radius>(particles[i]) = cell_radius;
        get<type>(particles[i]) = 0; // initially all cells are leaders

        //get<position>(p) = vdouble2(cell_radius,(i+1)*diameter); // x=2, uniformly in y
        get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=2, uniformly in y


    }

    particles.init_neighbour_search(vdouble2(0, 0), 5 * vdouble2(length_x, length_y), vbool2(false, false));


    /*
     * random initialisation
     */

//    for (int i=0; i<N; ++i) {
//        bool free_position = false;
//        particle_type::value_type p;
//        get<radius>(p) = cell_radius;
//        while(free_position == false){
//            get<position>(p) = vdouble2(cell_radius,uniform(gen)); // x=2, uniformly in y
//            free_position = true;
//            /*
//             * loop over all neighbouring particles within "diameter=2*radius" distance
//             */
//            for (auto tpl: euclidean_search(particles.get_query(),get<position>(p),diameter)) {
//                /*
//                 * tpl variable is a tuple containing:
//                 *  (0) -> neighbouring particle value_type
//                 *  (1) -> relative position of neighbouring particle
//                 *         from query point
//                 */
//                const vdouble2& dx = std::get<1>(tpl);
//                const particle_type::value_type& j = std::get<0>(tpl);
//                if (dx.norm() <  diameter) {
//                    free_position = false;
//                    break;
//                }
//            }
//        }
//        particles.push_back(p);
//    }

    // save particles before they move

    vtkWriteGrid("particles", t, particles.get_grid(true));


    // choose a set of random number between 0 and 2*pi, to avoid more rejections when it goes backwords (it would always be rejected)
    std::default_random_engine gen1;
    //gen1.seed(n_seed);
    std::uniform_real_distribution<double> uniformpi(0, 2 * M_PI);


    for (int t = 0; t < N_steps; t++) {



        // insert new cells at the start of the domain at insertion time (have to think about this insertion time)

        if (t % insertion_freq == 0) {
            bool free_position = false;
            particle_type::value_type f;
            get<radius>(f) = cell_radius;


            get<position>(f) = vdouble2(cell_radius, uniform(gen)); // x=2, uniformly in y
            free_position = true;
            /*
             * loop over all neighbouring leaders within "dem_diameter" distance
             */
            for (auto tpl = euclidean_search(particles.get_query(), get<position>(f), diameter); tpl != false; ++tpl) {

                vdouble2 diffx = tpl.dx();

                if (diffx.norm() < diameter) {
                    free_position = false;
                    break;
                }
            }

            get<type>(f) = 1; // new cells are followers
            //get<distance>(f) = 0;


            if (free_position) {
                particles.push_back(f);
                get<chain>(f) = 0;
                get<attached_to_id>(f) = 0;
            }

        }
        particles.update_positions();


        /////////////////////////////////////
        // grow domain


        if (t % freq_growth == 0) {

            //domain_length = domain_length + 1.0;

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
        //cout << "diff domain outside " << diff_domain << endl;




        // save the chemoattractant concentration with properly rescaled coordinates, UNIFORM DOMAIN GROWTH
//        for (int i = 0; i < length_x * length_y; i++) {
//            chemo_3col(i, 0) = chemo_3col_ind(i, 0) * (domain_length / length_x);
//        }
//        //cout << "domain length ratio " << domain_length/length_x << endl;

        // save the chemoattractant concentration with properly rescaled coordinates
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 0) = chemo_3col_ind(i, 0) * (domain_length / length_x);
        }


        // u column
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
        }


        // update chemoattractant profile

        /*
        * scaling factors for the same length scale
        * non uniform domain growth, only the first half grows
        * */

        double scaling_factor;

        // internalisation
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                //go through all the cells
                // leaders
                for (int k = 0; k < particles.size(); k++) {
                    vdouble2 x;
                    x = get<position>(particles[k]);

                    intern(i, j) = intern(i, j) + exp(-(((domain_length / length_x) * i - x[0]) *
                                                        ((domain_length / length_x) * i - x[0]) +
                                                        (j - x[1]) * (j - x[1])) /
                                                      (2 * R * R)); // mapping to fixed domain
                }
            }
        }

        for (int i = 1; i < length_x - 1; i++) {
            for (int j = 1; j < length_y - 1; j++) {


                // logistic production rate
                chemo_new(i, j) = dt * (D * ((1 / ((domain_length / length_x) * (domain_length / length_x))) *
                                             (chemo(i + 1, j) - 2 * chemo(i, j) + chemo(i - 1, j)) / (dx * dx) +
                                             (chemo(i, j + 1) - 2 * chemo(i, j) + chemo(i, j - 1)) / (dy * dy)) -
                                        (chemo(i, j) * lam / (2 * M_PI * R * R)) * intern(i, j) +
                                        kai * chemo(i, j) * (1 - chemo(i, j)) -
                                        double(domain_len_der) / double(domain_length) * chemo(i, j)) + chemo(i, j);

            }
        }

        // zero flux boundary conditions


        for (int i = 0; i < length_y; i++) {
            chemo_new(0, i) = chemo_new(1, i);
            chemo_new(length_x - 1, i) = chemo_new(length_x - 2, i);

        }

        for (int i = 0; i < length_x; i++) {
            chemo_new(i, 0) = chemo_new(i, 1);
            chemo_new(i, length_y - 1) = chemo_new(i, length_y - 2);
        }


        chemo = chemo_new; // update chemo concentration



        // save data to plot chemoattractant concentration
        ofstream output("matrix_growing_domain" + to_string(t) + ".csv");

        output << "x, y, z, u" << "\n" << endl;


        for (int i = 0; i < length_x * length_y; i++) {
            for (int j = 0; j < 4; j++) {
                output << chemo_3col(i, j) << ", ";
            }
            output << "\n" << endl;
        }


        /// update positions based on the domain growth, NON-uniform


        /// update positions uniformly based on the domain growth

        if (t % freq_growth == 0) {

            for (int i = 0; i < particles.size(); i++) {
                get<position>(particles)[i] *= vdouble2((domain_length / old_length), 1);
            }
            old_length = domain_length;
        }




        /////////////////////////////////////


        // Update positions based on the gradient



        // SEARCH MULTIPLE TIMES

        /*
         * create a random list of cell ids
         */


        int check_rep = 0; // check for repetitions, 0 no rep, 1 rep

        //for (int j=0; j<N_steps;j++){

        std::default_random_engine gen2;
        gen2.seed(t); // different seeds
        std::uniform_real_distribution<double> uniform_particles(0, particles.size()); // can only move forward

        VectorXi particle_id = VectorXi::Zero(particles.size());

        for (int i = 0; i < particles.size(); i++) {

            check_rep = 1; // set to 1 to enter the while loop
            while (check_rep == 1) {
                check_rep = 0; // it is initially zero and then will be changed to 1 if it is equivalent to others
                particle_id(i) = uniform_particles(gen2);


                for (int j = 0; j < i; j++) {
                    if (particle_id(i) == particle_id(j)) { check_rep = 1; }
                }
            }
            //cout << "ids " << particle_id(i) << endl;
        }


        // move all the leaders

        // pick a cell randomly

        for (int j = 0; j < particles.size(); j++) {


//            /*
// * phenotypic switching, based on chemoattractant concentration in front, +0.5
// * */
//
//            vdouble2 coord = get<position>(particles[particle_id(j)]);
//
//            // rescaled coord
//
//            double rescaled_coord;
//
//            rescaled_coord = (length_x / domain_length)*coord[0];
//
//
//
//            double chemo_in_front = chemo(round(rescaled_coord), round(coord[1]));
//            //cout << "chemo in front " << old_chemo << endl;
//
//
//            // if high concentration cells become leaders
//            if (chemo_in_front > chemo_leader ){
//
//                get<type>(particles[particle_id(j)]) = 0;
//            }
//            else{
//                get<type>(particles[particle_id(j)]) = 1;
//            }


            //go through all the leaders
            if (get<type>(particles[particle_id(j)]) == 0) {


                // go through cells in order
                //for (int i = 0; i < particles.size(); i++) {

                vdouble2 x;
                x = get<position>(particles[particle_id(j)]);


                /*
             * x_in variable will correspond to the coordinates on the non-updated domain (same as initial)
             * */


                double x_in; // x coordinate in initial domain length

                // Uniform domain growth


                x_in = (length_x / domain_length) * x[0];//
                l_filo_x = (length_x / domain_length) * l_filo_x; // rescale the length of filopodia as well


                if (t == 0) {
                    cout << "ratio" << length_x / domain_length << endl;
                }

                // create an array to store random directions
                array<double, filo_number +1> random_angle;
//            std::array<int, 3> sign_x;
//            std::array<int, 3> sign_y;
                for (int k = 0; k < filo_number+1; k++) {

                    double random_angle_tem = uniformpi(gen1);
//                int sign_x_tem, sign_y_tem;

                    // this is to have always inside the domain, but actually could sendi it out, and then the concentration there is zero
                    //while (round((x_in + sin(random_angle_tem) * l_filo_x)) < 0 ||
                    //       round((x_in + sin(random_angle_tem) * l_filo_x)) >
                    //      length_x - 1 || round(x[1] + cos(random_angle_tem) * l_filo_y) < 0 ||
                    //       round(x[1] + cos(random_angle_tem) * l_filo_y) > length_y - 1) {


                    random_angle_tem = uniformpi(gen1);


                    //}

                    random_angle[k] = random_angle_tem;

//                sign_x[j] = sign_x_tem;
//                sign_y[j] = sign_y_tem;


                }


                // choose which direction to move

                // if none of the criterion satisfied, do not move

                get<direction>(particles)[particle_id(j)] = vdouble2(0, 0);

                // store variables for concentration at new locations


                double old_chemo = chemo((round(x_in)), round(x)[1]);
                array<double, filo_number> new_chemo;


                for (int i = 0; i< filo_number; i++){

                    if (round((x_in + sin(random_angle[i]) * l_filo_x)) < 0 ||
                        round((x_in + sin(random_angle[i]) * l_filo_x)) >
                        length_x - 1 || round(x[1] + cos(random_angle[i]) * l_filo_y) < 0 ||
                        round(x[1] + cos(random_angle[i]) * l_filo_y) > length_y - 1) {
                        new_chemo[i] = 0;
                    } else {
                        new_chemo[i] = chemo(round((x_in + sin(random_angle[i]) * l_filo_x)),
                                             round(x[1] + cos(random_angle[i]) * l_filo_y));
                    }

                }


                // find maximum concentration of chemoattractant

                int chemo_max_number = 0;

                for (int i = 1; i < filo_number; i++){
                    if (new_chemo[chemo_max_number] < new_chemo[i]) {
                        chemo_max_number = i;
                    }
                }


                cout << "chemo max number " << chemo_max_number << endl;


                //if relative/absolute difference is smaller for the maximumu, move random direction
                //absolute
                //if (new_chemo[chemo_max_number] - old_chemo < diff_conc) {


                // relative
                if ((new_chemo[chemo_max_number] - old_chemo) / sqrt(old_chemo) < diff_conc) {

                    x += speed_l * vdouble2(sin(random_angle[3]), cos(random_angle[3]));
                    //cout << "print id " << id_[x] << endl;

                    x_in = (length_x / domain_length) * x[0];//uniform growth in the first part of the domain

                    //cout << "Position "<< x << endl;

                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move
                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                        }
                    }




                    // check that the position they want to move to is free and not out of bounds
                    if (free_position &&
                        (x_in) > 0 &&
                        (x_in) < length_x - 1 &&
                        (x[1]) > 0 &&
                        (x[1]) < length_y - 1) {
                        get<position>(particles)[particle_id(j)] += speed_l * vdouble2(sin(random_angle[filo_number]),
                                                                                       cos(random_angle[filo_number])); // update if nothing is in the next position
                        get<direction>(particles)[particle_id(j)] = speed_l * vdouble2(sin(random_angle[filo_number]),
                                                                                       cos(random_angle[filo_number]));
                    }

                } else {

                    x += speed_l * vdouble2(sin(random_angle[chemo_max_number]), cos(random_angle[chemo_max_number]));
                    //cout << "print id " << id_[x] << endl;

                    x_in = (length_x / domain_length) * x[0];//uniform growth in the first part of the domain

                    //cout << "Position "<< x << endl;

                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move
                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                        }
                    }


                    // check that the position they want to move to is free and not out of bounds
                    if (free_position &&
                        (x_in) > 0 &&
                        (x_in) < length_x - 1 &&
                        (x[1]) > 0 &&
                        (x[1]) < length_y - 1) {
                        get<position>(particles)[particle_id(j)] +=
                                speed_l * vdouble2(sin(random_angle[chemo_max_number]),
                                                   cos(random_angle[chemo_max_number])); // update if nothing is in the next position
                        get<direction>(particles)[particle_id(j)] =
                                speed_l * vdouble2(sin(random_angle[chemo_max_number]),
                                                   cos(random_angle[chemo_max_number]));
                    }


                }

            }

                //go through all the followers
                if (get<type>(particles[particle_id(j)]) == 1) {


                    vdouble2 x;
                    x = get<position>(particles[particle_id(j)]);

                    /*
                    * x_in variable will correspond to the coordinates on the non-updated domain (same as initial)
                    * */

                    double x_in; // x coordinate in initial domain length


                    // Non-uniform domain growth


                    x_in = (length_x / domain_length) * x[0];//uniform growth in the first part of the domain




                    // if the particle is part of the chain
                    if (get<chain>(particles[particle_id(j)]) > 0) {

                        get<direction>(particles)[particle_id(j)] = get<direction>(particles)[get<attached_to_id>(
                                particles[particle_id(j)])];

                        //try to move in the same direction as the cell it is attached to
                        vdouble2 x_chain = x + get<direction>(particles)[particle_id(j)];

                        // Non-uniform domain growth
                        double x_in_chain;

                        x_in_chain =
                                (length_x / domain_length) * x_chain[0];//uniform growth in the first part of the domain



                        bool free_position = true;

                        // update position of leaders, so that followers would have more space to move
                        //particles.update_positions();


                        for (auto pos = euclidean_search(particles.get_query(), x_chain, diameter);
                             pos != false; ++pos) {

                            //for (int i=0; i < particles.size(); i++) {
                            if (get<id>(*pos) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                //cout << "reject step " << get<id>(*pos) << endl;
                                free_position = false;
                                //get<chain>(particles)[particle_id(j)] = 0;
                            }
                        }


                        // check that the position they want to move to is free and not out of bounds, also it has to be non-zero, so that it would not be attached to non-moving cells
                        if (free_position &&
                            (x_in_chain) > 0 &&
                            (x_in_chain) < length_x - 1 && (x_chain[1]) > 0 &&
                            (x_chain[1]) < length_y - 1) {
                            //cout << "direction " << get<direction>(particles[particle_id(j)]) << endl;
                            get<position>(particles)[particle_id(j)] += get<direction>(particles[particle_id(j)]);

                        }


                        // check if it is not too far from the cell it was following

                        vdouble2 dist;

                        dist = get<position>(particles[particle_id(j)]) -
                               get<position>(particles[get<attached_to_id>(particles[particle_id(j)])]);

                        if (dist.norm() > l_filo_max) {
                            get<chain>(particles[particle_id(j)]) = 0;

                            // all the cells in the chain behind it become detached
                            for (int i = 0; i < particles.size(); ++i) {
                                if (get<chain_type>(particles[i]) == get<chain_type>(particles)[particle_id(j)]) {
                                    get<chain>(particles[i]) = 0;
                                }

                            }
                        }


                    }


                }



//go through all the followers
            if (get<type>(particles[particle_id(j)]) == 1){


                vdouble2 x;
                x = get<position>(particles[particle_id(j)]);

                /*
                * x_in variable will correspond to the coordinates on the non-updated domain (same as initial)
                * */

                double x_in; // x coordinate in initial domain length


                // Non-uniform domain growth


                x_in = (length_x / domain_length)*x[0];//uniform growth in the first part of the domain




                // if the particle is part of the chain
                if (get<chain>(particles[particle_id(j)]) > 0){

                    get<direction>(particles)[particle_id(j)] = get<direction>(particles)[get<attached_to_id>(particles[particle_id(j)])];

                    //try to move in the same direction as the cell it is attached to
                    vdouble2 x_chain = x + get<direction>(particles)[particle_id(j)];

                    // Non-uniform domain growth
                    double x_in_chain;

                    x_in_chain =
                            (length_x / domain_length) * x_chain[0];//uniform growth in the first part of the domain



                    bool free_position = true;

                    // update position of leaders, so that followers would have more space to move
                    //particles.update_positions();


                    for (auto pos = euclidean_search(particles.get_query(), x_chain, diameter); pos != false; ++pos) {

                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(*pos) !=
                            get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            //cout << "reject step " << get<id>(*pos) << endl;
                            free_position = false;
                            //get<chain>(particles)[particle_id(j)] = 0;
                        }
                    }


                    // check that the position they want to move to is free and not out of bounds, also it has to be non-zero, so that it would not be attached to non-moving cells
                    if (free_position &&
                        (x_in_chain) > 0 &&
                        (x_in_chain) < length_x - 1 && (x_chain[1]) > 0 &&
                        (x_chain[1]) < length_y - 1) {
                        //cout << "direction " << get<direction>(particles[particle_id(j)]) << endl;
                        get<position>(particles)[particle_id(j)] += get<direction>(particles[particle_id(j)]);

                    }


                    // check if it is not too far from the cell it was following

                    vdouble2 dist;

                    dist = get<position>(particles[particle_id(j)]) - get<position>(particles[get<attached_to_id>(particles[particle_id(j)])]);

                    if (dist.norm() > l_filo_max ){
                        get<chain>(particles[particle_id(j)]) = 0;

                        // all the cells in the chain behind it become detached
                        for (int i = 0 ; i< particles.size(); ++i){
                            if (get<chain_type>(particles[i]) == get<chain_type>(particles)[particle_id(j)]){
                                get<chain>(particles[i]) = 0;
                            }

                        }
                    }


                }



                // if the cell is not part of the chain
                    if (get<chain>(particles[particle_id(j)]) == 0) {


                        // check if there are any cells distance l_filo_y apart
                        // it can be either a leader or a follower already in a chain
                        //

                        // try to find a close leader
                        for (auto k = euclidean_search(particles.get_query(), x, l_filo_y); k != false; ++k) {


                            if (get<type>(*k) == 0) { // if it is close to a leader
                                //cout << "neighbours id leader" << get<id>(*k) << endl;
                                //get<direction>(particles[particle_id(j)]) = 0.2 * k.dx();
                                get<direction>(particles)[particle_id(j)] = get<direction>(*k);
                                get<chain>(particles)[particle_id(j)] = 1;
                                get<attached_to_id>(particles)[particle_id(j)] = get<id>(*k);
                                get<chain_type>(particles)[particle_id(j)] = get<id>(*k);
                            }

                        }


                        // try to find a close follower which is in a chain contact with a leader

                        for (auto k = euclidean_search(particles.get_query(), x, l_filo_y); k != false; ++k) {

                            //cout << "norm " << k.dx().norm() << endl;
                            //cout << "neighbours id " << get<id>(*k) << endl;

                            // if it is close to a follower that is part of the chain
                            if (get<type>(*k) == 1 && get<chain>(*k) > 0) {

                                if (get<id>(*k) != get<id>(particles[particle_id(j)])) {
                                    get<direction>(particles)[particle_id(j)] = get<direction>(*k);
                                    get<chain>(particles)[particle_id(j)] = get<chain>(*k) + 1;
                                    get<attached_to_id>(particles)[particle_id(j)] = get<id>(*k);
                                    get<chain_type>(particles)[particle_id(j)] = get<chain_type>(*k);
                                }


                            }

                        }






                        //try to move in the same direction as the cell it is attached to
                        vdouble2 x_chain = x + get<direction>(particles)[particle_id(j)];

                        // Non-uniform domain growth
                        double x_in_chain;

                        x_in_chain =
                                (length_x / domain_length) * x_chain[0];//uniform growth in the first part of the domain



                        bool free_position = true;

                        // update position of leaders, so that followers would have more space to move
                        //particles.update_positions();


                        for (auto pos = euclidean_search(particles.get_query(), x_chain, diameter);
                             pos != false; ++pos) {

                            //for (int i=0; i < particles.size(); i++) {
                            if (get<id>(*pos) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                //cout << "reject step " << get<id>(*pos) << endl;
                                free_position = false;
                                //get<chain>(particles)[particle_id(j)] = 0;
                            }
                        }




                        // check that the position they want to move to is free and not out of bounds, also it has to be non-zero, so that it would not be attached to non-moving cells
                        if (free_position &&
                            (x_in_chain) > 0 &&
                            (x_in_chain) < length_x - 1 && (x_chain[1]) > 0 &&
                            (x_chain[1]) < length_y - 1) {
                            //cout << "direction " << get<direction>(particles[particle_id(j)]) << endl;
                            get<position>(particles)[particle_id(j)] += get<direction>(particles[particle_id(j)]);

                        }


                        // if it hasn't found anything close, move randomly

                        if (get<chain>(particles[particle_id(j)]) == 0) {

                            double random_angle = uniformpi(gen1);

                            while (((x_in + sin(random_angle) * l_filo_y) < 0 ||
                                    ((x_in + sin(random_angle) * l_filo_y)) >
                                    length_x - 1 || (x[1] + cos(random_angle) * l_filo_y) < 0 ||
                                    (x[1] + cos(random_angle) * l_filo_y) > length_y - 1)) {
                                random_angle = uniformpi(gen1);
                            }


                            x += speed_f * vdouble2(sin(random_angle), cos(random_angle));
                            // Non-uniform domain growth


                            // Non-uniform domain growth, onl first half grows
                            // if in the first part of the domain
                            x_in = (length_x / domain_length) * x[0];//uniform growth in the first part of the domain


                            bool free_position = true; // check if the neighbouring position is free

                            // if this loop is entered, it means that there is another cell where I want to mov

                            for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {



                                //for (int i=0; i < particles.size(); i++) {
                                if (get<id>(*k) !=
                                    get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                    //cout << "reject step " << 1 << endl;
                                    free_position = false;
                                }
                            }
//                cout << " x coord " << round((x[0] * (length_x / domain_length))) << endl;
//                cout << "y coord " <<  round(x[1]) << endl;
//
//
//                if (round((x[0] * (length_x / domain_length))) > 0 &&
//                round((x[0] * (length_x / domain_length))) < length_x - 1 && round(x[1]) > 0 &&
//                round(x[1]) < length_y - 1){
//
//                    cout << "this condition satisfied" << endl;
//                }


                            // check that the position they want to move to is free and not out of bounds
                            if (free_position && x_in > 0 &&
                                (x_in) < length_x - 1 && (x[1]) > 0 &&
                                (x[1]) < length_y - 1) {
                                //cout << " moves " << endl;
                                //cout << "how frequently come in here " << endl;
                                get<position>(particles)[particle_id(j)] += speed_f * vdouble2(sin(random_angle),
                                                                                               cos(random_angle)); // update if nothing is in the next position
                            }

                        }

                    }


                    /*
                     * End of this new version
                     * */



//                for (auto k = euclidean_search(particles.get_query(), x, l_filo_y); k != false; ++k) {
//
//                    //cout << "norm " << k.dx().norm() << endl;
//                    //cout << "neighbours id " << get<id>(*k) << endl;
//
//                    // if it is close to a follower that is part of the chain
//                    if (get<type>(*k) == 1 && get<chain>(*k) == 1) {
//                        //cout << "neighbours id fol" << get<id>(*k) << endl;
//
//                        //get<direction>(particles[particle_id(j)]) = 0.1 * k.dx(); // move closer
//
//                        //get<position>(particles)[particle_id(j)] += speed_l * vdouble2(sin(random_angle[1]),
//                        //cos(random_angle[1])); // update if nothing is in the next position
//                        // check that it is not the same cell
//                        if (get<id>(*k) != get<id>(particles[particle_id(j)])){
//                            get<direction>(particles)[particle_id(j)] = get<direction>(*k);
//                            get<chain>(particles)[particle_id(j)] = 1;
//                            get<attached_to_id>(particles)[particle_id(j)] = get<id>(*k);
//                        }
//
//
//                    }
//
//                    if (get<type>(*k) == 0) { // if it is close to a leader
//                        //cout << "neighbours id leader" << get<id>(*k) << endl;
//                        //get<direction>(particles[particle_id(j)]) = 0.2 * k.dx();
//                        get<direction>(particles)[particle_id(j)] = get<direction>(*k);
//                        get<chain>(particles)[particle_id(j)] = 1;
//                        get<attached_to_id>(particles)[particle_id(j)] = get<id>(*k);
//                    }
//
//                }
//                // }
//
//
//                if (get<chain>(particles)[particle_id(j)] == 1) {
//                    get<direction>(particles)[particle_id(j)] = get<direction>(particles)[get<attached_to_id>(
//                            particles)[particle_id(j)]];
//
//
//                    // if part of the chain, same direction as the one that it follows, if not part of the chain, (0,0)
//                    vdouble2 x_chain = x + get<direction>(particles)[particle_id(j)];
//
//                    // Non-uniform domain growth
//                    double x_in_chain;
//
//                    x_in_chain =
//                            (length_x / domain_length) * x_chain[0];//uniform growth in the first part of the domain
//
//
//
//                    bool free_position = true;
//
//                    // update position of leaders, so that followers would have more space to move
//                    //particles.update_positions();
//
//
//                    for (auto pos = euclidean_search(particles.get_query(), x_chain, diameter); pos != false; ++pos) {
//
//                        //for (int i=0; i < particles.size(); i++) {
//                        if (get<id>(*pos) !=
//                            get<id>(particles[particle_id(j)])) { // check if it is not the same particle
//                            //cout << "reject step " << get<id>(*pos) << endl;
//                            free_position = false;
//                            //get<chain>(particles)[particle_id(j)] = 0;
//                        }
//                    }
//
//
//                    vdouble2 dir = get<direction>(particles[particle_id(j)]);
//
//
//
//                    // check that the position they want to move to is free and not out of bounds, also it has to be non-zero, so that it would not be attached to non-moving cells
//                    if (dir[0] != 0 && dir[1] != 0 && free_position && get<chain>(particles)[particle_id(j)] == 1 &&
//                        round(x_in_chain) > 0 &&
//                        round(x_in_chain) < length_x - 1 && round(x_chain[1]) > 0 &&
//                        round(x_chain[1]) < length_y - 1) {
//                        //cout << "direction " << get<direction>(particles[particle_id(j)]) << endl;
//                        get<position>(particles)[particle_id(j)] += get<direction>(particles[particle_id(j)]);
//
//                    }
//
//                    else{
//                        get<chain>(particles)[particle_id(j)] = 0; // it becomes dettached
//                        //get<attached_to_id>(particles)[particle_id(j)] = 0;
//                    }
//                }
//
//
//
//                if(get<chain>(particles)[particle_id(j)] == 0){
//
//                    double random_angle = uniformpi(gen1);
//
//
//                    while (round((x_in + sin(random_angle) * l_filo_x)) < 0 ||
//                           round((x_in + sin(random_angle)  * l_filo_x)) >
//                           length_x - 1 || round(x[1] + cos(random_angle) * l_filo_y) < 0 ||
//                           round(x[1] + cos(random_angle)  * l_filo_y) > length_y - 1) {
//
//                        random_angle = uniformpi(gen1);
//
//
//                    }
//
//                    x += speed_f * vdouble2(sin(random_angle), cos(random_angle));
//                    // Non-uniform domain growth
//
//
//                    // Non-uniform domain growth, onl first half grows
//                    // if in the first part of the domain
//                    x_in = (length_x / domain_length)*x[0];//uniform growth in the first part of the domain
//
//
//                    bool free_position = true; // check if the neighbouring position is free
//
//                    // if this loop is entered, it means that there is another cell where I want to mov
//
//                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {
//
//
//
//                        //for (int i=0; i < particles.size(); i++) {
//                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
//                            //cout << "reject step " << 1 << endl;
//                            free_position = false;
//                        }
//                    }
////                cout << " x coord " << round((x[0] * (length_x / domain_length))) << endl;
////                cout << "y coord " <<  round(x[1]) << endl;
////
////
////                if (round((x[0] * (length_x / domain_length))) > 0 &&
////                round((x[0] * (length_x / domain_length))) < length_x - 1 && round(x[1]) > 0 &&
////                round(x[1]) < length_y - 1){
////
////                    cout << "this condition satisfied" << endl;
////                }
//
//
//                    // check that the position they want to move to is free and not out of bounds
//                    if (free_position && x_in > 0 &&
//                        round(x_in) < length_x - 1 && round(x[1]) > 0 &&
//                        round(x[1]) < length_y - 1) {
//                        //cout << " moves " << endl;
//                        //cout << "how frequently come in here " << endl;
//                        get<position>(particles)[particle_id(j)] += speed_f * vdouble2(sin(random_angle), cos(random_angle)); // update if nothing is in the next position
//                    }
//                }
//
//                //get<chain>(particles[particle_id(j)]) == 0;
//                //cout << "direction of each particle " << get<direction>(particles[particle_id(j)]) << endl;
//                cout << "id of follower " << get<id>(particles[particle_id(j)]) << endl;

                    /*
                     * End of the old version
                     *
                     * */




                    /*
                     * Alternative phenotypic switching if a follower over takes a leader it becomes a leader and that leader follower.
                     * I will have to be careful when there will be channels because I will have to choose the closest leader
                     * */


                    // find the closest leader


                    // so that I would not go through all the cells I will choose the ones that are closer to the front

                    // minimum direction in x of the leaders

                    int min_index = 0;

                    for (int i = 1; i < N; ++i) {
                        if (get<position>(particles[i])[0] < get<position>(particles[min_index])[0]) {
                            min_index = i;
                        }

                    }


                    if (get<position>(particles[particle_id(j)])[0] > get<position>(particles[min_index])[0] + eps) {
                        // find distance to all the leaders
                        double distances[N];
                        vdouble2 dist_vector;
                        //check which one is the closest
                        for (int i = 0; i < N; ++i) {
                            // if it is one of the leaders

                            dist_vector = get<position>(particles[particle_id(j)]) - get<position>(particles[i]);
                            cout << "distance vector " << dist_vector << endl;
                            distances[i] = dist_vector.norm();
                            cout << "i = " << i << " distances " << distances[i] << endl;

                            for (int i = 0; i < N; ++i) {
                                cout << "i = " << i << " distances " << distances[i] << endl;
                            }

                            int winning_index = 0;
                            for (int i = 1; i < N; ++i) {
                                if (distances[i] < distances[winning_index]) {
                                    winning_index = i;
                                }
                            }

                            cout << "winning index " << winning_index << endl;

                            // if this closest leader is behind that follower, swap them
                            if (get<position>(particles[particle_id(j)])[0] >
                                get<position>(particles[winning_index])[0] + eps) {
                                particle_type::value_type tmp = particles[winning_index];


                                // their position swap

                                vdouble2 temp = get<position>(particles[winning_index]);
                                get<position>(particles[winning_index]) = get<position>(particles[particle_id(j)]);
                                get<position>(particles[particle_id(j)]) = temp;


                            }

                        }


                    }


                }
            }

            particles.update_positions();




            // save all time steps
#ifdef HAVE_VTK
            vtkWriteGrid("particles", t, particles.get_grid(true));
#endif


        }




//        /*
//         * return the density of cells in each of the fifth of the domain
//         */
//        const int domain_partition = int(domain_length / double(5));; // number of intervalas of 50 \mu m
//
//        VectorXi proportions = VectorXi::Zero(domain_partition); // integer with number of cells in particular part
//        //array<double, domain_partition> proportions;
//
//
//        double one_part = domain_length / double(domain_partition);
//
//        cout << "one part of the domain " << one_part << endl;
//
//        for (int i = 0; i < domain_partition; i++) {
//
//            cout << "number of cells " << particles.size() << endl;
//            for (int j = 0; j < particles.size(); j++) {
//                vdouble2 x = get<position>(particles[j]);
//                //cout<< "domain partition " << i*one_part << endl;
//                //cout << "x coordinate " << x[0] << endl;
//                cout << "position " << get<position>(particles[j]) << endl;
//                cout << "from here " << i * one_part << endl;
//                cout << "to here " << (i + 1) * one_part << endl;
//                if (i * one_part < x[0] && x[0] < (i + 1) * one_part) {
//                    proportions(i) += 1;
//                }
//            }
//
//        }
//
//
//        // for loops to count the number of cells in each of the fifth of the domain
//
//        return proportions;


    }



/*
 * main for futhest distance
 * */

//
//int main(){
//
//    const int number_parameters = 5; // parameter range
//
//    const int sim_num = 1;
//
//    MatrixXf all_distances = MatrixXf::Zero(number_parameters,sim_num); //matrix over which I am going to average
//
////n would correspond to different seeds
//    for (int n = 0; n < sim_num; n++) {
//
//        // define parameters that I will change
//        //VectorXf slope, threshold;
//        //array<double, 1> lam;
//
//        double speed_l[number_parameters];
//
//        for (int i = 0; i < number_parameters+1; ++i){
//            speed_l[i]= 0.1*(i+1);
//        }
//
//
//        VectorXf furthest_distance = VectorXf::Zero(number_parameters);
//
//#pragma omp parallel for
//        for (int i = 0; i < number_parameters; i++) {
//            for (int j = 0; j < 1; j++) {
//                cout << "iteration i " << i << endl;
//                furthest_distance(i, j) = distance_function(speed_l[i]);
//                cout << "parameters " << i << endl;
//            }
//            all_distances(i, n) = furthest_distance(i, 0);
//        }
//    }
////        /*
////        // save data to plot chemoattractant concentration
////        ofstream output("furthest_distance_matrix_cell_induced.csv");
////
////        MatrixXd furthest_distance_3col(number_parameters * number_parameters, 4), furthest_distance_3col_ind(number_parameters * number_parameters,
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
////        */
////
////
////
////
////        // This might be useful for matlab
////        ofstream output2("furthest_distance_matrix_matlab_cell_induced.csv");
////
////        for (int i = 0; i < number_parameters; i++) {
////
////            for (int j = 0; j < 1; j++) {
////
////                output2 << furthest_distance(i, j) << ", ";
////
////            }
////            output2 << "\n" << endl;
////        }
////
////    }
////
////
//    ofstream output3("simulations_cell_induced.csv");
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
//
//}





/*
 * main for proportions in different sections
 */


//// parameter analysis
//int main(){
//
//    const int number_parameters = 1; // parameter range
//    const int sim_num = 1;
//
//    VectorXi vector_check_length = proportions(0.005, 2); //just to know what the length is
//
//    int num_parts = vector_check_length.size(); // number of parts that I partition my domain
//
//    MatrixXf sum_of_all = MatrixXf::Zero(num_parts,number_parameters); // sum of the values over all simulations
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
//            threshold[i] = 0.5;
//            //threshold[i] = 0.005 * (i + 1);// 0.01;
//            //cout << "slope " << slope[i] << endl;
//
//        }
//
//
//
//
//        MatrixXi numbers = MatrixXi::Zero(num_parts,number_parameters); // can't initialise because do not know the size
//
//        cout << "stops here" << endl;
//
//#pragma omp parallel for
//        for (int i = 0; i < number_parameters; i++) {
//
//            //for (int j = 0; j < 1; j++) {
//
//            numbers.block(0,i,num_parts,1) = proportions(threshold[i], n);
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
//                sum_of_all(i,j) += numbers(i,j);
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