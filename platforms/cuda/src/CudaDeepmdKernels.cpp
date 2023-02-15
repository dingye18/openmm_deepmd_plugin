/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2018 Stanford University and the Authors.           *
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

#include "CudaDeepmdKernels.h"
#include "CudaDeepmdKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include <map>
#include <algorithm>

using namespace DeepmdPlugin;
using namespace OpenMM;
using namespace std;


CudaCalcDeepmdForceKernel::~CudaCalcDeepmdForceKernel(){
   return;
}

void CudaCalcDeepmdForceKernel::initialize(const System& system, const DeepmdForce& force){
    graph_file = force.getDeepmdGraphFile();
    type4EachParticle = force.getType4EachParticle();
    typesIndexMap = force.getTypesIndexMap();
    forceUnitCoeff = force.getForceUnitCoefficient();
    energyUnitCoeff = force.getEnergyUnitCoefficient();
    coordUnitCoeff = force.getCoordUnitCoefficient();

    lambda = force.getLambda();
    
    //natoms = system.getNumParticles();
    natoms = type4EachParticle.size();
    tot_atoms = system.getNumParticles();

    // Load the ordinary graph firstly.
    dp = DeepPot(graph_file);
    
    // Initialize the ordinary input and output array.
    // Initialize the input tensor.
    dener = 0.;
    dforce = vector<VALUETYPE>(natoms * 3, 0.);
    dvirial = vector<VALUETYPE>(9, 0.);
    dcoord = vector<VALUETYPE>(natoms * 3, 0.);
    dbox = vector<VALUETYPE>(9, 0.);
    //dtype = vector<int>(natoms, 0);    
    // Set atom type;
    //for(int ii = 0; ii < natoms; ii++){
        // ii is the atom index of each particle.
    //    dtype[ii] = typesIndexMap[type4EachParticle[ii]];
    //}

    for(std::map<int, string>::iterator it = type4EachParticle.begin(); it != type4EachParticle.end(); ++it){
        dp_particles.push_back(it->first);
        dtype.push_back(typesIndexMap[it->second]);
    }

    AddedForces = vector<double>(tot_atoms * 3, 0.0);
    // Set for CUDA context.
    cu.setAsCurrent();
    map<string, string> defines;
    defines["FORCES_TYPE"] = "double";
    networkForces.initialize(cu, 3*natoms, sizeof(double), "networkForces");
    CUmodule module = cu.createModule(CudaDeepmdKernelSources::DeepmdForce, defines);
    addForcesKernel = cu.getKernel(module, "addForces");
}


double CudaCalcDeepmdForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3> pos;
    context.getPositions(pos);
    Vec3 box[3];

    // Set box size.
    if (context.getSystem().usesPeriodicBoundaryConditions()){
        cu.getPeriodicBoxVectors(box[0], box[1], box[2]);
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
    // Assign the input coord for alchemical simulation.

    dp.compute (dener, dforce, dvirial, dcoord, dtype, dbox);

    // Transform the unit from output forces unit to KJ/(mol*nm)
    for(int ii = 0; ii < natoms; ii ++){
        int atom_index = dp_particles[ii];

        AddedForces[atom_index * 3 + 0] = lambda * dforce[ii * 3 + 0] * forceUnitCoeff;
        AddedForces[atom_index * 3 + 1] = lambda * dforce[ii * 3 + 1] * forceUnitCoeff;
        AddedForces[atom_index * 3 + 2] = lambda * dforce[ii * 3 + 2] * forceUnitCoeff;
    }
    dener = dener * energyUnitCoeff;

    if (includeForces) {
        // Change to OpenMM CUDA context.
        cu.setAsCurrent();
        networkForces.upload(AddedForces);
        int paddedNumAtoms = cu.getPaddedNumAtoms();
        void* args[] = {&networkForces.getDevicePointer(), &cu.getForce().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer(), &natoms, &paddedNumAtoms};
        cu.executeKernel(addForcesKernel, args, natoms);
    }
    return dener;
}



