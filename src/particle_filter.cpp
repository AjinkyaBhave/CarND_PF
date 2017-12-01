/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std_pos[]) {
	// Number of particles in filter
	n_particles = 100;
	// Random number generator for Gaussian distribution
	default_random_engine rand_gen;
	// Gaussian distribution for x, y and theta
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	// Temporary class to store particle initial data
	Particle p;
	// Initialize all particles to first position (based on estimates of 
	// x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// Particle positions are initialised in map coordinate frame.
	for (int i=0; i<n_particles; i++){
		p.x = x+dist_x(rand_gen);
		p.y = y+dist_y(rand_gen);
		p.theta = theta+dist_theta(rand_gen);
		p.id = i;
		p.weight = 1;
		particles.push_back(p);
	}
	// Filter state is now initialised
	is_initialized = true;
	cout << "Filter initialised with " << particles.size() << " particles.\n"; 
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// Random number generator for Gaussian distribution
	default_random_engine rand_gen;
	// Gaussian distribution for x, y and theta
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	
	for (int i=0; i<particles.size(); i++){
		// Predict next state using bicycle model 
		particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta));
		particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t));
		particles[i].theta = particles[i].theta + yaw_rate*delta_t;
		
		// Add Gaussian noise to state
		particles[i].x = particles[i].x + dist_x(rand_gen);
		particles[i].y = particles[i].y + dist_y(rand_gen);
		particles[i].theta = particles[i].theta + dist_theta(rand_gen);	
	}
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	std::vector<LandmarkObs> predicted;
	LandmarkObs tempObs;
	double min_obs_error;
	double obs_error;
	int id;
	
	for(int i=0; i<particles.size(); i++){
		predicted.clear();
		for(int j=0; j<observations.size();j++){
			tempObs.x = particles[i].x + cos(particles[i].theta)*observations[i].x - sin(particles[i].theta)*observations[i].y;
			tempObs.y = particles[i].y + sin(particles[i].theta)*observations[i].x + cos(particles[i].theta)*observations[i].y;
			min_obs_error = sensor_range;
			for(int k=0; k<map_landmarks.landmark_list.size(); k++){
				obs_error = (map_landmarks.landmark_list[i].x_f-tempObs.x)*(map_landmarks.landmark_list[i].x_f-tempObs.x);
				obs_error += (map_landmarks.landmark_list[i].y_f-tempObs.y)*(map_landmarks.landmark_list[i].y_f-tempObs.y);
				obs_error = sqrt(obs_error);
				if(obs_error < min_obs_error){
					min_obs_error = obs_error;
					id = map_landmarks.landmark_list[i].id_i;
				}
			}
			tempObs.id = id;
			predicted.push_back(tempObs);
			
		}
		
	}
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    // particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
