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

#include "ReferenceDeepmdKernels.h"
#include "DeepmdForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include <typeinfo>
#include <iostream>
#include <map>
#include <algorithm>
#include <limits>

using namespace DeepmdPlugin;
using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static Vec3* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (Vec3*) data->periodicBoxVectors;
}

ReferenceCalcDeepmdForceKernel::~ReferenceCalcDeepmdForceKernel(){return;}

void ReferenceCalcDeepmdForceKernel::initialize(const System& system, const DeepmdForce& force) {
    graph_file = force.getDeepmdGraphFile();
    type4EachParticle = force.getType4EachParticle();
    typesIndexMap = force.getTypesIndexMap();
    forceUnitCoeff = force.getForceUnitCoefficient();
    energyUnitCoeff = force.getEnergyUnitCoefficient();
    coordUnitCoeff = force.getCoordUnitCoefficient();
    lambda = force.getLambda();
    natoms = type4EachParticle.size();
    tot_atoms = system.getNumParticles();

    // Fetch parameters for adaptive DP region.
    isFixedRegion = force.isFixedRegion();
    if (!isFixedRegion){
        center_atoms = force.getCenterAtoms();
        radius = force.getRegionRadius();
        atom_names4dp_forces = force.getAtomNames4DPForces();
        sel_num4type = force.getSelNum4EachType();
        topology = force.getTopology();
    }
    
    // Set up the dtype and input, output arrays.
    if(isFixedRegion){
        dener = 0.;
        dforce = vector<VALUETYPE>(natoms * 3, 0.);
        dvirial = vector<VALUETYPE>(9, 0.);
        dcoord = vector<VALUETYPE>(natoms * 3, 0.);
        dbox = vector<VALUETYPE>(9, 0.);

        for(std::map<int, string>::iterator it = type4EachParticle.begin(); it != type4EachParticle.end(); ++it){
            dp_particles.push_back(it->first);
            dtype.push_back(typesIndexMap[it->second]);
        }    
    } else {
        natoms = 0;
        for(map<string, int>::iterator it = sel_num4type.begin(); it != sel_num4type.end(); ++it){
            cum_sum4type[it->first] = vector<int>(2, 0);
            cum_sum4type[it->first][0] = natoms;
            natoms += it->second;
            cum_sum4type[it->first][1] = natoms;
            for (int i = 0; i < it->second; i++){
                dtype.push_back(typesIndexMap[it->first]);
            }
        }

        dener = 0.;
        dforce = vector<VALUETYPE>(natoms * 3, 0.);
        dvirial = vector<VALUETYPE>(9, 0.);
        dcoord = vector<VALUETYPE>(natoms * 3, 0.);
        dbox = {}; // Empty vector for adaptive region.
        daparam = vector<int>(natoms, 0.);
    }

    // Initialize DeepPot.
    dp = DeepPot(graph_file);
    AddedForces = vector<double>(tot_atoms * 3, 0.0);

}

double ReferenceCalcDeepmdForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& force = extractForces(context);

    if (isFixedRegion){
        // Set box size.
        if (context.getSystem().usesPeriodicBoundaryConditions()){
            Vec3* box = extractBoxVectors(context);
            // Transform unit from nanometers to angstrom.
            dbox[0] = box[0][0] * coordUnitCoeff;
            dbox[1] = box[0][1] * coordUnitCoeff;
            dbox[2] = box[0][2] * coordUnitCoeff;
            dbox[3] = box[1][0] * coordUnitCoeff;
            dbox[4] = box[1][1] * coordUnitCoeff;
            dbox[5] = box[1][2] * coordUnitCoeff;
            dbox[6] = box[2][0] * coordUnitCoeff;
            dbox[7] = box[2][1] * coordUnitCoeff;
            dbox[8] = box[2][2] * coordUnitCoeff;
        }else{
            dbox = {}; // No PBC.
        }
        // Set input coord.
        for(int ii = 0; ii < natoms; ++ii){
            // Multiply by coordUnitCoeff means the transformation of the unit from nanometers to required input unit for positions in trained DP model.
            int atom_index = dp_particles[ii];
            dcoord[ii * 3 + 0] = pos[atom_index][0] * coordUnitCoeff;
            dcoord[ii * 3 + 1] = pos[atom_index][1] * coordUnitCoeff;
            dcoord[ii * 3 + 2] = pos[atom_index][2] * coordUnitCoeff;
        }
        dp.compute (dener, dforce, dvirial, dcoord, dtype, dbox);

    } else {
        
        vector<int> addOrNot; // Whether to add the dp forces for selected atoms. 
        map<string, vector<int>> dp_region_atoms = DeepmdPlugin::SearchAtomsInRegion
        (pos, center_atoms, radius, topology, atom_names4dp_forces, addOrNot);
        
        for(map<string, vector<int>>::iterator it = dp_region_atoms.begin(); it != dp_region_atoms.end(); ++it){
            string atom_type = it->first;
            vector<int> atom_indices = it->second;
            int atom_num = atom_indices.size();
            if (atom_num > sel_num4type[atom_type]){
                throw OpenMMException("The number of atoms in the adaptive region is larger than the number of atoms selected for DP forces.");
            }
            int cum_sum_start = cum_sum4type[atom_type][0];
            int cum_sum_end = cum_sum4type[atom_type][1];
            for (int i = 0; i < atom_num; i++){
                int atom_index = atom_indices[i];
                int dp_index = cum_sum_start + i;
                dcoord[dp_index * 3 + 0] = pos[atom_index][0] * coordUnitCoeff;
                dcoord[dp_index * 3 + 1] = pos[atom_index][1] * coordUnitCoeff;
                dcoord[dp_index * 3 + 2] = pos[atom_index][2] * coordUnitCoeff;
                daparam[dp_index] = 1;
            }
        }
        vector<VALUETYPE> dfparam = {};
        
        dp.compute(dener, dforce, dvirial, dcoord, dtype, dfparam, daparam);

        for(int ii = 0; ii < addOrNot.size(); ii++){
            if (addOrNot[ii] == 0){
                dforce[ii * 3 + 0] = 0.;
                dforce[ii * 3 + 1] = 0.;
                dforce[ii * 3 + 2] = 0.;
            }
        }
    }
    
    

    // Transform the unit from output forces unit to KJ/(mol*nm)
    for(int ii = 0; ii < natoms; ii ++){
        int atom_index = dp_particles[ii];
        AddedForces[atom_index * 3 + 0] = lambda * dforce[ii * 3 + 0] * forceUnitCoeff;
        AddedForces[atom_index * 3 + 1] = lambda * dforce[ii * 3 + 1] * forceUnitCoeff;
        AddedForces[atom_index * 3 + 2] = lambda * dforce[ii * 3 + 2] * forceUnitCoeff;
    }
    dener = lambda * dener * energyUnitCoeff;

    if(includeForces){
        for(int ii = 0; ii < tot_atoms; ii ++){
        force[ii][0] += AddedForces[ii * 3 + 0];N
        force[ii][1] += AddedForces[ii * 3 + 1];
        force[ii][2] += AddedForces[ii * 3 + 2];
        }
    }
    // Return energy.
    return dener;
}


