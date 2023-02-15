/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "DeepmdForce.h"
#include "internal/DeepmdForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <sys/stat.h>

using namespace DeepmdPlugin;
using namespace OpenMM;
using namespace std;

inline bool exists(const std::string& name) {
    struct stat buffer;   
    return (stat (name.c_str(), &buffer) == 0); 
}


DeepmdForce::DeepmdForce(const string& GraphFile){
    graph_file  = GraphFile;
    if (!exists(graph_file)){
        throw OpenMMException("Graph file not found: "+graph_file);
    }
    // Initialize dp model
    DeepPot tmp_dp = DeepPot(graph_file);
    this->numb_types = tmp_dp.numb_types();
    this->cutoff = tmp_dp.cutoff();
    tmp_dp.get_type_map(this->type_map);
}

DeepmdForce::DeepmdForce(const string& GraphFile, const double& lambda){
    graph_file  = GraphFile;
    this->lambda = lambda;
    if (!exists(graph_file)){
        throw OpenMMException("Graph file not found: "+graph_file);
    }
    // Initialize dp model
    DeepPot tmp_dp = DeepPot(graph_file);
    this->numb_types = tmp_dp.numb_types();
    this->cutoff = tmp_dp.cutoff();
    tmp_dp.get_type_map(this->type_map);
}

DeepmdForce::~DeepmdForce(){
    type4EachParticle.clear();
    particleGroup4EachType.clear();
    typesIndexMap.clear();
}


void DeepmdForce::setPBC(const bool use_PBC){
    // By default, use_pbc is set to be true.
    use_pbc = use_PBC;
}

void DeepmdForce::setUnitTransformCoefficients(const double coordCoefficient, const double forceCoefficient, const double energyCoefficient){
    coordCoeff = coordCoefficient;
    forceCoeff = forceCoefficient;
    energyCoeff = energyCoefficient;
}

double DeepmdForce::getCoordUnitCoefficient() const {return coordCoeff;}
double DeepmdForce::getForceUnitCoefficient() const {return forceCoeff;}
double DeepmdForce::getEnergyUnitCoefficient() const {return energyCoeff;}

double DeepmdForce::getCutoff() const {return cutoff;}
int DeepmdForce::getNumberTypes() const {return numb_types;}
string DeepmdForce::getTypesMap() const {return type_map;}

const string& DeepmdForce::getDeepmdGraphFile() const{return graph_file;}
const map<int, string>& DeepmdForce::getType4EachParticle() const{return type4EachParticle;}
const map<string, vector<int>>& DeepmdForce::getParticles4EachType() const{return particleGroup4EachType;}
const map<string, int>& DeepmdForce::getTypesIndexMap() const{return typesIndexMap;}

void DeepmdForce::addParticle(const int particleIndex, const string particleType){
    auto insertResult = type4EachParticle.insert(pair<int, string>(particleIndex, particleType));
    if(insertResult.second == false){
        throw OpenMMException("Failed to add particle, duplicate key.");
    }
    auto it = particleGroup4EachType.find(particleType);
    if (it == particleGroup4EachType.end()){
        particleGroup4EachType[particleType] = vector<int>();
        particleGroup4EachType[particleType].push_back(particleIndex);
    }else{
        particleGroup4EachType[particleType].push_back(particleIndex);
    }
}

void DeepmdForce::addType(const int typeIndex, const string Type){
    auto it = typesIndexMap.find(Type);
    if(it == typesIndexMap.end()){
        typesIndexMap[Type] = typeIndex;
    }else{
        if(typeIndex != it->second){
            throw OpenMMException("Type Index duplicated.");
        }
    }
}

void DeepmdForce::addBond(const int particle1, const int particle2){
    bondsList.push_back(make_pair(particle1, particle2));
}

const vector<pair<int, int>> DeepmdForce::getBondsList() const{
    return bondsList;
}


ForceImpl* DeepmdForce::createImpl() const {
    return new DeepmdForceImpl(*this);
}

void DeepmdForce::updateParametersInContext(Context& context) {
    // Nothing to be done here.
    return;
}

void DeepmdForce::setLambda(const double lambda){
    this->lambda = lambda;
}

double DeepmdForce::getLambda() const {return lambda;}
