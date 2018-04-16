
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

int main() {


    // model parameters

    int length_x = 30;//240;//40;//4; // length of the chemoattractant vector, for fixed domain
    double domain_length = 30; //this variable is for the actual domain length
    double old_length = 30;// this is important for the update of the positions
    const int length_y = 12;//120;//20;//4;
    double cell_radius = 0.75;//0.5; // radius of a cell
    const double diameter = 2 * cell_radius;//2 // diameter in which there have to be no cells, equivalent to size of the cell
    const int N_steps = 800; // number of times the cells move up the gradient
    const size_t N = 4; // initial number of cells
    double l_filo = 27.5/10;//2; // sensing radius
    double diff_conc = 0.05; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 1; // determines how frequently domain grows (actually not relevant because it will go every timestep)
    int insertion_freq = 1;
    double speed_l = 0.5;//0.05; // speed of a leader cell
    double speed_f = 0.5;//0.08; // speed of a follower cell
    double dettach_prob = 0.5; // probability that a follower cell which is on trail looses the trail
    double chemo_leader = 0.9; //0.5; // phenotypic switching happens when the concentration of chemoattractant is higher than this (presentation video 0.95), no phenotypic switching
    double eps = 1; // for phenotypic switching, the distance has to be that much higher
    double track_spacing = speed_f *2; // spacing between positions on the track
    int track_length = 50;

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
    double a = 0.008;//0.008;//0.23/10;
    double L_inf = 86.76;
    double t_s = 16;//4.31*10;
    double constant = 29.12;

    double domain_len_der = 0; // initialise derivative of the domain growth function

    // parameters for the dynamics of chemoattractant concentration


    double D = 1; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0; // initialise time, redundant
    double dt = 0.000001; // time step
    double dx = 1; // space step in x direction, double to be consistent with other types
    double dy = 1; // space step in y direction
    double kai = 1/100;//0.0001/10; // to 1 /h production rate of chemoattractant


    // parameters for internalisation

    double R = cell_radius;//7.5/10; // \nu m cell radius
    double lam = 500;//(100)/10; // to 1000 /h chemoattractant internalisation


    // matrix that stores the values of concentration of chemoattractant
    MatrixXf chemo = MatrixXf::Zero(length_x, length_y);
    MatrixXf chemo_new = MatrixXf::Zero(length_x, length_y);

    // initialise internalisation matrix
    MatrixXf intern = MatrixXf::Zero(length_x, length_y);

    for (int i = 0;i<length_x;i++){
        for (int j = 0; j< length_y; j++){
            chemo(i,j) = 1; // uniform concentration initially
            chemo_new(i,j) = 1; // this is for later updates
        }    // generate initial chemoattractant concentration

    }


    // four columns for x, y, z, u (z is necessary for paraview)

    // form a matrix which would store x,y,z,u

    MatrixXf chemo_3col(length_x*length_y,4), chemo_3col_ind(length_x*length_y,2); // need for because that is how paraview accepts data, third dimension is just zeros

    // x, y coord, 1st and 2nd columns respectively
    int k = 0;
    // it has to be 3D for paraview
    while (k<length_x*length_y){
        for (int i = 0;i<length_x;i++){
            for (int j = 0; j< length_y; j++){
                chemo_3col_ind(k,0) = i;
                chemo_3col_ind(k,1) = j;
                chemo_3col(k,2) = 0;
                k += 1;
            }
        }
    }


    // y and x (initially) column
    for (int i=0;i<length_x*length_y;i++){
        chemo_3col(i,1) = chemo_3col_ind(i,1);
        chemo_3col(i,0) = chemo_3col_ind(i,0);
    }



    // create an array to keep track of all the positions of leader cells

    vdouble2 track_position [track_length][N] = {vdouble2(0,0)};
    for (int i =0; i < track_length; i++){
        for(int j = 0; j < N; j++){
            track_position[i][j] = vdouble2(0,0);
        }
    }
    int track_time[N] = {0}; // vector that stores the time values for each leader when there was a sufficiently big change in the position


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
    ABORIA_VARIABLE(attached_leader_nr,int,"attached_leader_nr")
    ABORIA_VARIABLE(attached_at_time_step,int,"attached_at_time_step")
    ABORIA_VARIABLE(type,int,"type") // 0 if a cell is a leader, 1 if follower


    ABORIA_VARIABLE(chain, double, "chain") // stores whether attached to a leader or follower
    // stores the distance to the closest neighbour, if less than thresold
    typedef Particles<std::tuple<radius, type, attached_to_id, attached_leader_nr, attached_at_time_step, direction, chain>, 2> particle_type; // 2 stands for dimension


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
    std::uniform_real_distribution<double> uniform(2,length_y-1);

    /*
     * initialise neighbour search with 2d cuboid domain,
     * periodic in x and y
     */

    particles.init_neighbour_search(vdouble2(0,0), 5*vdouble2(length_x,length_y), vbool2(false,false));


    /*
     * compact initialisation
     */

    for (int i=0; i<N; ++i) {


        get<radius>(particles[i]) = cell_radius;
        get<type>(particles[i]) = 0; // initially all cells are leaders

        //get<position>(p) = vdouble2(cell_radius,(i+1)*diameter); // x=2, uniformly in y
        get<position>(particles[i]) = vdouble2(cell_radius,(i+1)*double(length_y-1)/double(N)-0.5 * double(length_y-1)/double(N)); // x=2, uniformly in y


    }

    particles.init_neighbour_search(vdouble2(0,0), 5*vdouble2(length_x,length_y), vbool2(false,false));


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

    vtkWriteGrid("particles",t,particles.get_grid(true));


    // choose a set of random number between 0 and 2*pi, to avoid more rejections when it goes backwords (it would always be rejected)
    std::default_random_engine gen1;
    std::uniform_real_distribution<double> uniformpi(0,2*M_PI);


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

                if (diffx.norm() <  diameter) {
                    free_position = false;
                    break;
                }
            }

            get<type>(f) = 1; // new cells are followers
            get<attached_at_time_step>(f) == -1;
            get<attached_leader_nr>(f) == -1;


            if (free_position ) {
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

                // x_in = (length_x / domain_length)*x[0];

                // Non-uniform domain growth


                x_in = (length_x / domain_length)*x[0];//uniform growth in the first part of the domain




                // create an array to store random directions
                std::array<double, 3> random_angle;
//            std::array<int, 3> sign_x;
//            std::array<int, 3> sign_y;
                for (int k = 0; k < 3; k++) {

                    double random_angle_tem = uniformpi(gen1);
//                int sign_x_tem, sign_y_tem;

                    while (round((x_in + sin(random_angle_tem) * l_filo)) < 0 ||
                           round((x_in + sin(random_angle_tem) * l_filo)) >
                           length_x - 1 || round(x[1] + cos(random_angle_tem) * l_filo) < 0 ||
                           round(x[1] + cos(random_angle_tem) * l_filo) > length_y - 1) {
                        random_angle_tem = uniformpi(gen1);

//                    if (sin(random_angle_tem) < 0) {
//                        sign_x_tem = -1;
//                    } else { sign_x_tem = 1; }
//
//                    if (cos(random_angle_tem) < 0) {
//                        sign_y_tem = -1;
//                    } else { sign_y_tem = 1; }

                    }

                    random_angle[k] = random_angle_tem;

//                sign_x[j] = sign_x_tem;
//                sign_y[j] = sign_y_tem;


                }


                // choose which direction to move

                // if none of the criterion satisfied, do not move

                get<direction>(particles)[particle_id(j)] = vdouble2(0, 0);

                // store variables for concentration at new locations


                double old_chemo = chemo((round(x_in)), round(x)[1]);

                double new_chemo_1 = chemo(round((x_in + sin(random_angle[0]) * l_filo)),
                                           round(x[1] + cos(random_angle[0]) * l_filo));

                double new_chemo_2 = chemo(round((x_in + sin(random_angle[1]) * l_filo)),
                                           round(x[1] + cos(random_angle[1]) * l_filo));


                //if both smaller, move random direction
                //absolute
                if (new_chemo_1 - old_chemo < diff_conc && new_chemo_2 - old_chemo < diff_conc) {


                    // relative
                    //if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) < diff_conc && (new_chemo_2- old_chemo)/sqrt(old_chemo) < diff_conc){

                    x += speed_l * vdouble2(sin(random_angle[2]), cos(random_angle[2]));
                    //cout << "print id " << id_[x] << endl;

                    x_in = (length_x / domain_length)*x[0];//uniform growth in the first part of the domain

                    //cout << "Position "<< x << endl;

                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move
                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                            break;
                        }
                    }




                    // check that the position they want to move to is free and not out of bounds
                    if (free_position &&
                        round(x_in) >= 0 &&
                        round(x_in) <= length_x - 1 &&
                        round(x[1]) >= 0 &&
                        round(x[1] ) <= length_y - 1) {
                        get<position>(particles)[particle_id(j)] += speed_l * vdouble2(sin(random_angle[2]),
                                                                                       cos(random_angle[2])); // update if nothing is in the next position
                        get<direction>(particles)[particle_id(j)] = speed_l * vdouble2(sin(random_angle[2]),
                                                                                       cos(random_angle[2]));
                    }

                }
                    //cout << "stops here " << endl;
                    // if first direction greater, second smaller
                    //absolute
                else if (new_chemo_1 - old_chemo > diff_conc && new_chemo_2 - old_chemo < diff_conc) {

                    //relative
                    //else if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) > diff_conc && (new_chemo_2 - old_chemo)/sqrt(old_chemo) < diff_conc){

                    x += speed_l * vdouble2(sin(random_angle[0]), cos(random_angle[0]));
                    //cout << "print id " << id_[x] << endl;

                    // Non-uniform domain growth


                    x_in = (length_x / domain_length)*x[0];//uniform growth in the first part of the domain


                    //cout << "Position "<< x << endl;
                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move

                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                            break;
                        }
                    }


                    // check that the position they want to move to is free and not out of bounds
                    if (free_position &&
                        round(x_in) >= 0 &&
                        round(x_in) <= length_x - 1 &&
                        round(x[1] ) >= 0 &&
                        round(x[1]) <= length_y - 1) {
                        get<position>(particles)[particle_id(j)] += speed_l * vdouble2(sin(random_angle[0]),
                                                                                       cos(random_angle[0])); // update if nothing is in the next position
                        get<direction>(particles)[particle_id(j)] = speed_l * vdouble2(sin(random_angle[0]),
                                                                                       cos(random_angle[0]));
                    }

                }
                    // if first smaller, second bigger

                    //absolute
                else if (new_chemo_1 - old_chemo < diff_conc && new_chemo_2 - old_chemo > diff_conc) {

                    //relative
                    //else if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) < diff_conc && (new_chemo_2 - old_chemo)/sqrt(old_chemo) > diff_conc){



                    x += speed_l * vdouble2(sin(random_angle[1]), cos(random_angle[1]));
                    //cout << "print id " << id_[x] << endl;

                    // Non-uniform domain growth


                    x_in = (length_x / domain_length)*x[0];//uniform growth in the first part of the domain


                    //cout << "Position "<< x << endl;
                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move
                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                            break;
                        }
                    }



                    // check that the position they want to move to is free and not out of bounds
                    if (free_position &&
                        round(x_in) >= 0 &&
                        round(x_in) <= length_x - 1 &&
                        round(x[1]) >= 0 &&
                        round(x[1]) <= length_y - 1) {
                        get<position>(particles)[particle_id(j)] += speed_l * vdouble2(sin(random_angle[1]),
                                                                                       cos(random_angle[1])); // update if nothing is in the next position
                        get<direction>(particles)[particle_id(j)] = speed_l * vdouble2(sin(random_angle[1]),
                                                                                       cos(random_angle[1]));
                    }
                    //break;
                }
                    // if both greater choose the bigger one

                    // absolute
                else if (new_chemo_1 - old_chemo > diff_conc && new_chemo_2 - old_chemo > diff_conc) {

                    //relative
                    //else if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) > diff_conc && (new_chemo_2 - old_chemo)/sqrt(old_chemo) > diff_conc){


                    // if first is greater than the second
                    if (new_chemo_1 > new_chemo_2) {
                        x += speed_l * vdouble2(sin(random_angle[0]), cos(random_angle[0]));
                        //cout << "print id " << id_[x] << endl;
                        // Non-uniform domain growth


                        x_in = (length_x / domain_length)*x[0];//uniform growth in the first part of the domain


                        //cout << "Position "<< x << endl;
                        bool free_position = true; // check if the neighbouring position is free

                        // if this loop is entered, it means that there is another cell where I want to move
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            //for (int i=0; i < particles.size(); i++) {
                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                //cout << "reject step " << 1 << endl;
                                free_position = false;
                                break;
                            }
                        }



                        // check that the position they want to move to is free and not out of bounds
                        if (free_position &&
                            round(x_in) >= 0 &&
                            round(x_in) <= length_x - 1 &&
                            round(x[1] ) >= 0 &&
                            round(x[1]) <= length_y - 1) {
                            get<position>(particles)[particle_id(j)] += speed_l * vdouble2(sin(random_angle[0]),
                                                                                           cos(random_angle[0])); // update if nothing is in the next position
                            get<direction>(particles)[particle_id(j)] = speed_l * vdouble2(sin(random_angle[0]),
                                                                                           cos(random_angle[0]));

                        }

                    }
                        // if second is greater than the first
                    else if (new_chemo_1 < new_chemo_2) {
                        x += speed_l * vdouble2(sin(random_angle[1]), cos(random_angle[1]));
                        //cout << "print id " << id_[x] << endl;
                        // Non-uniform domain growth


                        x_in = (length_x / domain_length)*x[0];//uniform growth in the first part of the domain


                        //cout << "Position "<< x << endl;
                        bool free_position = true; // check if the neighbouring position is free

                        // if this loop is entered, it means that there is another cell where I want to move
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            //for (int i=0; i < particles.size(); i++) {
                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                //cout << "reject step " << 1 << endl;
                                free_position = false;
                                break;
                            }
                        }


                        // check that the position they want to move to is free and not out of bounds
                        if (free_position &&
                            round(x_in) >= 0 &&
                            round(x_in) <= length_x - 1 &&
                            round(x[1]) >= 0 &&
                            round(x[1]) <= length_y - 1) {
                            get<position>(particles)[particle_id(j)] += speed_l * vdouble2(sin(random_angle[1]),
                                                                                           cos(random_angle[1])); // update if nothing is in the next position
                            get<direction>(particles)[particle_id(j)] = speed_l * vdouble2(sin(random_angle[1]),
                                                                                           cos(random_angle[1]));

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



                get<chain>(particles[particle_id(j)]) = 0;
                get<attached_to_id>(particles[particle_id(j)]) = -1;

                // initially set direction to zero, if not attached to a leader
//                if (get<chain>(particles)[particle_id(j)] == 0){
//                    get<direction>(particles[particle_id(j)]) = vdouble2(0,0);
//                }

                //update direction, if a particle is part of a chain
//                if(get<chain>(particles)[particle_id(j)] == 1) {
//                    get<direction>(particles[particle_id(j)]) = get<direction>(particles[get<attached_to_id>(particles[particle_id(j)])]);
//                }

                //bool free_position = true; // check if the neighbouring position is free
                // if a particle is not part of the chain, look for cells that are leaders or part of a chain
                //if (get<chain>(particles[particle_id(j)]) == 0) {


                cout << "current id " << get<id>(particles)[particle_id(j)] << endl;

                for (auto k = euclidean_search(particles.get_query(), x, l_filo); k != false; ++k) {

                    //cout << "norm " << k.dx().norm() << endl;
                    //cout << "neighbours id " << get<id>(*k) << endl;

                    // if it is close to a follower that is part of the chain
                    if (get<type>(*k) == 1 && get<chain>(*k) == 1) {
                        //cout << "neighbours id fol" << get<id>(*k) << endl;

                        //get<direction>(particles[particle_id(j)]) = 0.1 * k.dx(); // move closer

                        //get<position>(particles)[particle_id(j)] += speed_l * vdouble2(sin(random_angle[1]),
                        //cos(random_angle[1])); // update if nothing is in the next position
                        // check that it is not the same cell
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])){
                            get<direction>(particles)[particle_id(j)] = get<direction>(*k);
                            get<chain>(particles)[particle_id(j)] = 1;
                            get<attached_to_id>(particles)[particle_id(j)] = get<id>(*k);
                        }


                    }

                    if (get<type>(*k) == 0) { // if it is close to a leader
                        //cout << "neighbours id leader" << get<id>(*k) << endl;
                        //get<direction>(particles[particle_id(j)]) = 0.2 * k.dx();
                        get<direction>(particles)[particle_id(j)] = get<direction>(*k);
                        get<chain>(particles)[particle_id(j)] = 1;
                        get<attached_to_id>(particles)[particle_id(j)] = get<id>(*k);
                    }

                }
                // }

                // if part of the chain, same direction as the one that it follows, if not part of the chain, (0,0)
                vdouble2 x_chain = x + get<direction>(particles)[particle_id(j)];

                // Non-uniform domain growth
                double x_in_chain;

                x_in_chain = (length_x / domain_length)*x_chain[0];//uniform growth in the first part of the domain



                bool free_position = true;

                // update position of leaders, so that followers would have more space to move
                //particles.update_positions();


                for (auto pos = euclidean_search(particles.get_query(), x_chain, diameter); pos != false; ++pos) {

                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(*pos) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                        //cout << "reject step " << get<id>(*pos) << endl;
                        free_position = false;
                        break;
                        //get<chain>(particles)[particle_id(j)] = 0;
                    }
                }


                vdouble2 dir = get<direction>(particles[particle_id(j)]);



                // check that the position they want to move to is free and not out of bounds, also it has to be non-zero, so that it would not be attached to non-moving cells
                if (dir[0] != 0 && dir[1] != 0 && free_position && get<chain>(particles)[particle_id(j)] == 1 && round(x_in_chain) > 0 &&
                    round(x_in_chain) < length_x - 1 && round(x_chain[1]) > 0 &&
                    round(x_chain[1]) < length_y - 1 ) {
                    //cout << "direction " << get<direction>(particles[particle_id(j)]) << endl;
                    get<position>(particles)[particle_id(j)] += get<direction>(particles[particle_id(j)]);

                }



                else{
                    get<chain>(particles)[particle_id(j)] = 0; // it becomes dettached
                    //get<attached_to_id>(particles)[particle_id(j)] = 0;
                }

                if(get<chain>(particles)[particle_id(j)] == 0) {

                    int found_track = 0;

                    /*
                     * Check if there is a track nearby
                     * */

                    // check what the closest neighbour is

                    vdouble2 diff; // difference in position
                    double previous_diff = dist_thres; // initialise distance threshold

                    /*for all the particles that are still not part of the chain
                    try to find a track of a leader cell to which it is close to
                     but make sure I do not
                     */
                    for (int k = 0; k < N; k++) {
                        for (int i = 0; i < track_length; i++) {
                            diff = x - track_position[i][k];
                            if (diff.norm() < dist_thres) {
                                //cout << "previous diff " << previous_diff << endl;
                                if (diff.norm() < previous_diff) {
                                    for (int f = 0;
                                         f <
                                         particles.size(); f++) { // make sure other particles are not in that position

                                        //cout << " time step value " << get<attached_at_time_step>(particles[f]) << endl;
                                        //cout << " attached to a leader " << get<attached_leader_nr>(particles[f])

                                        //cout << "value of j " << j << endl;
                                        //cout << "value of k " << k << endl;
                                        if (get<attached_at_time_step>(particles[f]) != j &&
                                            get<attached_leader_nr>(particles[f]) != k) {
                                            // have to make sure that it is not the same position as the leader
                                            //for (int p = 0; p < particles.size(); p++) {

                                            // this was in order to avoid follower cells moving on top of leader cells
                                            /* if ((track_position[j][k] - get<position>(particles[p])).norm() < 1) {
                                                 get<chain>(particles[i]) = 0;
                                                 break;
                                             }*/

                                            get<attached_at_time_step>(particles[particle_id(j)]) = i;
                                            get<attached_leader_nr>(particles[particle_id(j)]) = k;
                                            //}
                                        }
                                    }
                                }
                            }//else { get<chain>(particles[i]) = 0; }
                        }
                    }


                    vdouble2 x_can = track_position[get<attached_at_time_step>(
                            particles[particle_id(j)])+1][get<attached_leader_nr>(
                            particles[particle_id(j)])]; //position where the cell wants to move

                    // check if that position is free
                    //cout << "Position "<< x << endl;

                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move
                    for (auto k = euclidean_search(particles.get_query(), x_can, diameter); k != false; ++k) {


                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                            break;
                            cout << "how frequently free position is false" << endl;
                        }
                    }



                    // it is possible that randomly they get dettached


                    // if the cell is part of the chain update its position
                    if (free_position &&
                        round((x_can[0] * (length_x / domain_length))) >= 0 &&
                        round((x_can[0] * (length_x / domain_length))) < length_x - 1 && round(x_can[1]) > 0 &&
                        round(x_can[1]) < length_y - 1) {

                        get<attached_at_time_step>(particles[particle_id(j)]) += 1;

                        get<position>(particles[particle_id(j)]) = track_position[get<attached_at_time_step>(
                                particles[particle_id(j)])][get<attached_leader_nr>(particles[particle_id(j)])];
                        found_track = 1;
                        cout << "found track" << endl;
                    }

                    // if it hasn't found a track
                    if (found_track == 0) {


                        double random_angle = uniformpi(gen1);


                        while (round((x_in + sin(random_angle) * l_filo)) < 0 ||
                               round((x_in + sin(random_angle) * l_filo)) >
                               length_x - 1 || round(x[1] + cos(random_angle) * l_filo) < 0 ||
                               round(x[1] + cos(random_angle) * l_filo) > length_y - 1) {

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
                                break;
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
                        if (free_position && round(x_in) >= 0 &&
                            round(x_in) < length_x - 1 && round(x[1]) > 0 &&
                            round(x[1]) < length_y - 1) {
                            //cout << " moves " << endl;
                            //cout << "how frequently come in here " << endl;
                            get<position>(particles)[particle_id(j)] += speed_f * vdouble2(sin(random_angle),
                                                                                           cos(random_angle)); // update if nothing is in the next position
                        }
                    }

                    //get<chain>(particles[particle_id(j)]) == 0;
                    //cout << "direction of each particle " << get<direction>(particles[particle_id(j)]) << endl;
                    cout << "id of follower " << get<id>(particles[particle_id(j)]) << endl;

                    /*
                     * Alternative phenotypic switching if a follower over takes a leader it becomes a leader and that leader follower.
                     * I will have to be careful when there will be channels because I will have to choose the closest leader
                     * */


                                    // find the closest leader

                for (auto k = euclidean_search(particles.get_query(), get<position>(particles[particle_id(j)]), 2*diameter); k != false; ++k) {

                    if (get<type>(*k) == 0 && get<position>(particles[particle_id(j)])[0] > get<position>(*k)[0] + eps){
                        get<type>(particles[particle_id(j)]) = 0; // that particle becomes a leader
                        get<type>(*k) = 1; // leader becomes a follower
                        // their ids swap as well
                        int temp_id = get<id>(*k);
                        get<id>(*k) = get<id>(particles[particle_id(j)]);
                        get<id>(particles[particle_id(j)]) = temp_id;
                    }
                }



                    // in order find the closest leader
                    int clos_number = 0;

//                    for (auto k = euclidean_search(particles.get_query(), get<position>(particles[particle_id(j)]),
//                                                   2 * diameter); k != false; ++k) {
//
//                        int temp_id;
//
//                        if (get<type>(*k) == 0 &&
//                            get<position>(particles[particle_id(j)])[0] > get<position>(*k)[0] + eps) {
//                            // their ids swap as well
//                            temp_id = get<id>(*k); // temporarily id of sufficiently close leader
//                            clos_number = 1;
//                        }
//
//                        // if I have already found a close leader, check if it is closer to that follower
//
//
//                        if (clos_number == 1){
//
//                            vdouble2 dist_vector;
//                            vdouble2 dist_vector_previous;
//                            double dist;
//
//                            dist_vector = get<position>(particles[particle_id(j)]) - get<position>(*k);
//                            dist_vector_previous = get<position>(particles[particle_id(j)]) - get<position>(particles[temp_id]);
//
//                            if (dist_vector.norm()< dist_vector_previous.norm()){
//                             temp_id = get<id>(*k);
//                            }
//                            clos_number = 0;
//
//                        }
//
//                        if (clos_number == 0){
//                                get<type>(particles[particle_id(j)]) = 0; // that particle becomes a leader
//                                get<type>(particles[temp_id]) = 1; // leader becomes a follower
//                                // their ids swap as well
//                                get<id>(particles[temp_id]) = get<id>(particles[particle_id(j)]);
//                                get<id>(particles[particle_id(j)]) = temp_id;
//                                clos_number = 1;
//                        }
//
//                    }


                }

            }


        }

        particles.update_positions();

        // store the last track_length positions of all the leaders if the distance is sufficiently big

        cout << "track position " << track_position[0][1] << endl;




            //cout << 'track time ' << track_time[i] << endl;
            // check if new position is sufficiently far, but also not too far
            for(int i=0; i<particles.size(); i++) {
                if (get<type>(particles[i]) == 0){
                    vdouble2 diff;
                    int id_lead; // need this because leaders are first four ids
                    id_lead = get<id>(particles[i]);
                    // if shorter than track length
                    if (track_time[id_lead] < track_length) {
                        diff = track_position[track_time[id_lead]][id_lead] - get<position>(particles[id_lead]); // the difference is some intermediate vecotr value
                        cout << "track position " << track_position[track_time[id_lead]][id_lead] << endl;

                    }// if longer than track length
                    else {
                        diff = track_position[track_length-1][id_lead] - get<position>(particles[id_lead]); // the last position
                        cout << "track position last " << track_position[track_length][id_lead] << endl;
                    }
                    cout << "get position " << get<position>(particles[id_lead]) << endl;
                    cout << "norm difference " << diff.norm() << endl;
                    if (diff.norm() > track_spacing) {
                        track_time[id_lead] += 1; // update time steps

                        if (track_time[id_lead] < track_length) {
                            cout << "is it ever her " << endl;
                            track_position[track_time[id_lead]][id_lead] = get<position>(particles[id_lead]);
                        }// if longer than track length
                        if (track_time[id_lead] > track_length) {
                            for (int j = 0; j < track_length; j++) {
                                track_position[j][id_lead] = track_position[j + 1][id_lead];// vector shifts by one
                            }
                            track_position[track_length-1][id_lead] = get<position>(particles[id_lead]); // new position added

                        }
                    }
                    cout << "track time " << track_time[id_lead]<< endl;
                }

            }





        // save all time steps
#ifdef HAVE_VTK
        vtkWriteGrid("particles", t, particles.get_grid(true));
#endif


    }

    // save data to track positions
    ofstream output1("track_positions.csv");

    output1 << "x, y, z, u" << "\n" << endl;


    for (int i=0;i<track_length;i++){
        for(int j=0;j<2;j++){
            if(j==1){
                output1 << 0 << ", ";
            }
            else{output1 << track_position[i][j] << ", ";}
        }
        output1 << "\n" << endl;
    }




}
