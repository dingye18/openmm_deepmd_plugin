/* -------------------------------------------------------------------------- *
 *                                OpenMM-Deepmd                               *
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
#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace DeepmdPlugin;
using namespace OpenMM;
using namespace deepmd;
using namespace std;

extern "C" void registerDeepmdSerializationProxies();


void testSerialization() {
    const double TOL = 1e-5;
    const string graph = "../tests/frozen_model/water.pb";
    const double coordUnitCoeff = 10;
    const double forceUnitCoeff = 964.8792534459;
    const double energyUnitCoeff = 96.48792534459;
    const double temperature = 300;

    // Create a Force.
    DeepmdForce dp_force = DeepmdForce(graph);
    
    stringstream buffer;
    XmlSerializer::serialize<DeepmdForce>(&dp_force, "Force", buffer);
    DeepmdForce* copy = XmlSerializer::deserialize<DeepmdForce>(buffer);

    // Compare the two forces to see if they are identical.
    DeepmdForce& dp_force2 = *copy;
    ASSERT_EQUAL(dp_force2.getDeepmdGraphFile(), dp_force.getDeepmdGraphFile());
    return;
}

int main() {
    try {
        registerDeepmdSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
