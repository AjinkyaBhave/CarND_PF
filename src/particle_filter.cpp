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
	// Define maximum size for particles and weights vectors
	particles.reserve(n_particles);
	weights.reserve(n_particles);
	// Random number generator for Gaussian distribution
	std::random_device rand_dev;
	std::mt19937 rand_gen(rand_dev());
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
		weights.push_back(p.weight);
	}
	// Filter state is now initialised
	is_initialized = true;
	cout << "Filter initialised with " << particles.size() << " particles.\n"; 
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// Random number generator for Gaussian distribution
	std::random_device rand_dev;
	std::mt19937 rand_gen(rand_dev());
	// Gaussian distribution for x, y and theta
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	
	for (int i=0; i<particles.size(); i++){
		// Predict next state using bicycle model with large yaw rate
		if(fabs(yaw_rate) > 0.0001){
			particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t));
			particles[i].theta = particles[i].theta + yaw_rate*delta_t;	
		}
		// Vehicle is going in a straight line
		else{
			particles[i].x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			particles[i].y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
		}
		
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
	// Update the weights of each particle using a multi-variate Gaussian distribution. 
	// The observations are given in the VEHICLE'S coordinate system. Your particles are located
	// according to the MAP'S coordinate system. You will need to transform between the two systems.
	
	// Store detected landmarks within sensor range for a particle
	std::vector<LandmarkObs> predicted;
	// Temporary landmark data storage
	LandmarkObs temp_obs;
	// Temporary variables for distance error between observation and landmark
	double min_obs_error = sensor_range;
	double obs_error = sensor_range;
	int id = 0;
	// Multi-variate gaussian probability of all observations given particle pose
	double posterior_prob = 1;
	// Sum of all weights for normalisation
	//double sum_weights = 0;
	
	// Associate each particle with correct observed landmarks
	for(int i=0; i<particles.size(); i++)
	{
		// Empty list for current particle
		predicted.clear();
		// Initialise for weight calculation later
		posterior_prob = 1;
		for(int j=0; j<observations.size();j++)
		{
			// Convert from vehicle to world coordinate frame
			temp_obs.x = particles[i].x + cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y;
			temp_obs.y = particles[i].y + sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y;
			min_obs_error = sensor_range;
			for(int k=0; k<map_landmarks.landmark_list.size(); k++)
			{
				if(dist(map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, particles[i].x, particles[i].y) <= sensor_range){
					obs_error = dist(map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, temp_obs.x, temp_obs.y); 
					if(obs_error < min_obs_error)
					{
						min_obs_error = obs_error;
						id = map_landmarks.landmark_list[k].id_i;
					} //if()
				} // if()				
			} // for()
			temp_obs.id = id;
			predicted.push_back(temp_obs);
			posterior_prob *= gaussian(temp_obs, map_landmarks.landmark_list[id-1], std_landmark);
		} // for()
		// Assign posterior probability of all observations as weight of particle
		particles[i].weight = posterior_prob;
		weights[i] = particles[i].weight;		
		//sum_weights = sum_weights + particles[i].weight;
	} //for()

		/*// Normalise particle weights to create valid probability distribution
	for (int i=0;i<particles.size();i++){
		particles[i].weight = particles[i].weight/sum_weights;
		weights[i] = particles[i].weight;
	}*/
	
}

void ParticleFilter::resample() {
	//Resample particles with replacement with probability proportional to their weights. 
	
	// Random number generator for discrete weighted distribution
	std::random_device rand_dev;
	std::mt19937 rand_gen(rand_dev());
	// Create a discrete distribution with particle weights
	std::discrete_distribution<> discrete_prob(weights.begin(), weights.end());
	// Vector of new particles drawn with replacement based on weights
	std::vector<Particle> new_particles(n_particles);
	// Index of sampled particle
	int idx;
	
	for (int i=0;i<n_particles;i++)
	{
		idx = discrete_prob(rand_gen);
		//new_particles.push_back(particles[idx]);
		new_particles[i] = particles[idx];
	}
	particles = new_particles;
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
