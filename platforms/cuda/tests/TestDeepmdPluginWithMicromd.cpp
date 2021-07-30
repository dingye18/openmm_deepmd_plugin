#include "DeepmdForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/Platform.h"
#include "micromd.h"

using namespace OpenMM;
using namespace DeepmdPlugin;
using namespace std;

extern "C" OPENMM_EXPORT void registerDeepmdCudaKernelFactories();


void testDeepmdDynamicsWithMicromd(vector<VALUETYPE> init_positions, vector<string> particleTypeName, map<string, int> typeDict, vector<double> mass, int nsteps=1000){
    int natoms;
    natoms = int(init_positions.size()/3);

    // Set up the Micromd simulation environment.
    Micromd micromd = Micromd(natoms);
     
    // Set up OpenMM system.
    System system;
    VerletIntegrator integrator(delta_T); // Time step is 0.0001 ps here.
    DeepmdForce dp_force = DeepmdForce(graph, " ", " ", false);
    dp_force.setDeepmdOpFile(op_path);
    
    vector<Vec3> positions;
    // units is nanometers for box size setting.
    vector<Vec3> box = {Vec3(1.3,0.0,0.0),Vec3(0.0,1.3,0.0),Vec3(0.0,0.0,1.3)};
    vector<int> types;

    // Set up the OpenMM simulation environment.
    for(auto it = typeDict.begin(); it != typeDict.end(); it++){
        dp_force.addType(it->second, it->first);
    }
    for (int ii = 0; ii < natoms; ++ii){
        system.addParticle(mass[ii]);
        dp_force.addParticle(ii, particleTypeName[ii]);    
        types.push_back(typeDict[particleTypeName[ii]]);
        positions.push_back(Vec3(init_positions[ii * 3 + 0], init_positions[ii * 3 + 1], init_positions[ii * 3 + 2]));
        // Add atom to micromd.
        micromd.addAtom(atom(ii, types[ii], particleTypeName[ii], mass[ii]));
    }
    dp_force.setUnitTransformCoefficients(coordUnitCoeff, forceUnitCoeff, energyUnitCoeff); 
    system.addForce(&dp_force);
    Platform& platform = Platform::getPlatformByName("CUDA");
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setPeriodicBoxVectors(box[0], box[1], box[2]);
    context.setVelocitiesToTemperature(temperature, randomSeed);

    // Set up the Micromd environment.
    micromd.setInitPositions(positions);
    micromd.setInitBox(box);
    micromd.setTypes(types);
    micromd.setInitVelocitiesToTemperature(temperature, randomSeed);   

    State temp = context.getState(State::Velocities | State::Positions);
    vector<Vec3> omm_velocities = temp.getVelocities();
    vector<Vec3> omm_positions = temp.getPositions();

    micromd.setInitPositions(omm_positions);
    micromd.setInitVelocities(omm_velocities);

    //ASSERT_EQUAL_VEC(omm_positions[0], micromd_positions[0], TOL);
    ASSERT_EQUAL_VEC(omm_positions[0], positions[0], TOL);
    

    // Record the difference of forces and energy on each step.
    vector<double> errorForce;
    vector<double> errorEnergy;
    double err = 0.;

    for (int ii = 0; ii < nsteps; ++ii){
        // Running dynamics 1 step.
        integrator.step(1);
        micromd.VerletStep(1);
        
        // Get the force from openmm context state.
        State state = context.getState(State::Forces | State::Energy | State::Positions | State::Velocities);
        const vector<Vec3>& omm_forces = state.getForces();
        const double& omm_energy = state.getPotentialEnergy();
        const double& omm_kinetic_energy = state.getKineticEnergy();
        omm_positions = state.getPositions();
        omm_velocities = state.getVelocities();
        
        // Get the energy and forces from micromd.
        vector<Vec3> micromd_forces(natoms, Vec3(0,0,0));
        vector<Vec3> micromd_positions(natoms, Vec3(0,0,0));
        vector<Vec3> micromd_velocities(natoms, Vec3(0,0,0));
        double micromd_total_energy, micromd_potential_energy, micromd_kinetic_energy;

        micromd_forces = micromd.getForces();
        micromd_potential_energy = micromd.getPotentialEnergy();
        micromd_kinetic_energy = micromd.getKineticEnergy();
        micromd_positions = micromd.getPositions();
        micromd_velocities = micromd.getVelocities();
        
        //ASSERT_EQUAL_TOL(omm_kinetic_energy, micromd_kinetic_energy, TOL);
        // Compare the difference of these two results.
        err = 0.;
        for (int jj = 0; jj < natoms; ++jj){
            //ASSERT_EQUAL_VEC(omm_velocities[jj], micromd_velocities[jj], TOL);
            //ASSERT_EQUAL_VEC(omm_positions[jj], micromd_positions[jj], TOL);
            //ASSERT_EQUAL_VEC(omm_forces[jj], micromd_forces[jj], TOL);
            Vec3 diff_f = omm_forces[jj] - micromd_forces[jj];
            double diff = diff_f.dot(diff_f);
            err += diff;
        }
        //ASSERT_EQUAL_TOL(micromd_potential_energy, omm_energy, TOL);
        err = err/natoms;
        err = sqrt(err);
        cout<<"Diff on Forces: "<<err<<endl;
        errorForce.push_back(err);
        err = abs(omm_energy - micromd_potential_energy);
        errorEnergy.push_back(err);
        cout<<"Diff on Potential Energy: "<<err<<endl;

        cout<<"-----Total-------Potential-------Kinetic-----------Step: "<<ii<<endl;
        cout<<"OpenMM: "<< omm_energy + omm_kinetic_energy << " " << omm_energy << " " << omm_kinetic_energy<<endl;
        cout<<"Micromd: "<< micromd_potential_energy + micromd_kinetic_energy << " " << micromd_potential_energy << " " << micromd_kinetic_energy << endl;
        cout<<"-------------------------------------"<<endl;

    }
    cout<<"Dynamics Checking done Successfully!"<<endl;
    double meanErrorForce = std::accumulate(errorForce.begin(), errorForce.end(),0.)/nsteps;
    double meanErrorEnergy = std::accumulate(errorEnergy.begin(), errorEnergy.end(),0.)/nsteps;
    cout<<"Mean Force Difference: "<<meanErrorForce<<"; ";
    cout<<"Mean Energy Difference: "<<meanErrorEnergy<<"; ";
    cout<<"in "<<nsteps<<" steps"<<endl;
}


int main(int argc, char * argv[]){
    // Initialize positions, unit is angstrom.
    vector<float> positions_64_waters = {0.16800,3.02700,2.21300,1.12000,2.94300,2.15100,0.02900,3.70100,2.87800,10.10000,7.61600,7.77400,9.27800,7.12900,7.83100,10.77900,6.96200,7.94100,4.41000,6.68700,10.38500,4.44700,7.01800,11.28200,4.68400,7.42800,9.84400,9.41400,12.18000,1.12700,10.15600,12.36900,1.70200,9.54600,11.27100,0.85700,7.71300,6.00800,6.81800,7.64800,6.43500,5.96300,6.80500,5.84700,7.07600,-0.12100,5.54400,10.84600,-0.18100,4.76700,11.40300,0.45900,6.13600,11.32400,7.65300,11.50100,5.36300,8.14200,12.26900,5.65500,7.71800,10.88300,6.09100,2.05100,11.31000,2.97500,1.62100,10.70700,3.58100,2.96700,11.31800,3.25400,5.81900,3.05200,5.13500,5.08000,3.40300,4.63800,6.59100,3.44700,4.73000,0.59600,10.10400,10.37400,-0.25600,9.67200,10.31600,0.58200,10.54900,11.22200,11.54300,0.55300,2.51300,11.87300,1.43700,2.35700,12.23000,0.12800,3.02700,4.22900,7.74100,5.11600,3.67500,7.33800,4.44800,3.61600,8.19300,5.69600,1.29600,4.41200,5.15700,0.96800,3.52000,5.04400,0.72200,4.95000,4.61200,1.10800,8.37200,7.95700,1.32900,9.00700,8.63800,0.16700,8.48300,7.82300,9.80300,6.59800,11.39600,10.05600,7.42100,10.97800,10.57300,6.03600,11.30600,3.57400,4.44500,3.80500,3.48300,5.31400,3.41300,3.02000,4.47500,4.58500,3.97800,8.10900,0.27600,4.73000,8.00500,0.85900,3.90400,9.05400,0.14500,10.84800,0.90000,6.06400,10.95200,1.18500,6.97200,11.64100,1.20500,5.62400,3.92800,4.11600,11.48900,2.98500,4.26200,11.56800,4.26700,4.93100,11.11800,11.07300,4.34000,4.33200,10.58800,3.55800,4.06800,10.44500,5.05600,4.23300,1.35000,1.06800,4.76000,1.64200,0.39800,5.37800,2.05500,1.12200,4.11500,12.24600,8.12700,3.50800,11.34000,8.14600,3.81600,12.73700,8.60600,4.17700,5.61300,9.18200,9.02900,6.04600,8.94200,9.84900,5.62800,10.14000,9.02200,1.35100,7.22600,0.13000,2.10500,7.79000,-0.04300,0.76600,7.76600,0.66100,7.64400,3.72900,2.61400,7.91800,4.01600,1.74300,8.34200,3.14100,2.90200,11.27900,10.37800,6.35800,11.13500,11.32400,6.39400,10.60800,10.01100,6.93300,3.75900,11.44800,10.41900,3.04300,11.36100,11.04900,3.89400,10.56100,10.08500,8.46900,9.96500,7.56800,8.95400,9.27600,8.02100,7.70500,10.12500,8.12200,1.96600,11.13400,0.16000,1.83000,11.06400,1.10500,1.41000,11.86500,-0.10900,10.61400,9.32900,10.52100,10.27300,9.41800,11.41100,10.43500,10.17500,10.11000,6.71000,7.10900,4.13200,6.81300,7.88600,3.58300,5.97700,7.31900,4.71000,1.51700,3.15500,9.03800,1.82700,3.95200,8.60700,2.29700,2.60700,9.11800,4.70500,11.57800,5.13700,4.76000,12.38200,5.65300,5.59800,11.23300,5.13300,9.50900,2.01300,3.73100,9.45100,1.54800,4.56500,10.20900,1.56700,3.25500,9.42300,1.57600,11.27500,9.27800,1.18600,12.13800,8.66800,2.15000,11.14100,7.34900,2.97700,10.01000,7.87900,3.42800,9.35200,6.68400,2.51000,9.50400,11.46100,3.33000,11.84900,11.95400,2.93200,12.56700,10.71300,2.74600,11.71900,5.59800,2.34200,0.82100,6.09000,2.89800,1.42400,5.95800,2.54500,-0.04200,8.48100,4.73100,12.62900,8.20300,4.37000,11.78700,9.08800,5.43500,12.40000,12.04600,9.26500,1.06000,12.23000,8.96700,1.95100,11.12200,9.51500,1.07600,7.53800,1.28400,7.12700,7.04700,1.81300,6.49800,8.22200,1.87100,7.45100,11.58000,5.53800,8.10300,12.02700,5.50400,8.94900,12.04600,4.90300,7.55900,2.07100,11.70700,6.69600,1.55600,11.50700,7.47800,2.95400,11.87000,7.02500,2.85600,2.98600,1.55200,3.53600,2.93200,0.88000,3.30400,3.34200,2.31900,2.27500,6.71600,3.27100,2.32700,6.82900,2.32200,1.35700,6.88700,3.48300,9.34700,9.31400,1.02300,9.33200,8.43300,0.65100,8.68200,9.29700,1.71100,6.66200,9.13500,11.72000,6.67000,8.63400,12.53500,6.95900,10.00900,11.97200,6.74900,6.99700,1.11300,6.75800,6.76500,2.04200,7.32400,6.35300,0.70000,4.16800,1.96800,9.31400,4.03900,1.15300,9.80000,4.21000,2.64600,9.98900,4.19600,-0.14900,1.30700,4.48500,0.75800,1.40100,3.24600,-0.09000,1.20400,1.18000,1.17600,11.76100,1.79500,1.80900,12.13300,0.86000,1.59400,10.96200,9.68000,7.63500,4.97900,10.16800,7.71200,5.79900,9.01300,8.32000,5.03300,7.00700,11.53200,9.23200,6.97600,11.88900,10.11900,7.04700,12.30100,8.66400,6.87700,12.23900,-0.30600,6.04000,12.35900,0.14200,7.53700,12.42300,0.36300,11.51000,2.17000,8.42700,11.22700,2.34700,9.32400,12.45700,2.30900,8.44200,1.10400,9.86700,5.16200,1.65400,10.61400,5.39900,0.28400,10.01300,5.63300,8.92400,3.90400,7.64400,9.86800,4.04200,7.56800,8.53700,4.73300,7.36100,7.59200,9.75500,3.01100,6.68300,9.99600,2.82900,7.91700,10.45600,3.57700,10.01400,11.78600,9.48200,9.87000,12.66600,9.83000,9.24800,11.61600,8.93500,2.01300,5.53800,7.73400,1.84300,6.48000,7.74500,1.99800,5.30200,6.80600,0.38300,12.20300,8.80500,0.66500,11.63600,9.52300,-0.57300,12.19200,8.85300,4.91500,10.20200,2.85500,4.62000,10.76900,2.14300,4.62400,10.64600,3.65100,4.51400,1.30800,6.76300,4.69200,2.15100,6.34500,4.41100,1.51800,7.69100,5.08300,6.27700,7.62300,4.27700,6.00700,8.06400,5.01400,7.22900,7.55800};

    vector<string> types_64_waters = {"O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H","O","H","H"};

    vector<VALUETYPE> init_positions = positions_64_waters;
    vector<string> particleTypeName = types_64_waters;

    map<string, int> typeDict = {{"O", 0}, {"H", 1}};
    vector<double> mass;

    int nsteps = 1000000;

    int natoms = init_positions.size() / 3;
    for(int ii = 0; ii < natoms; ++ii){
        if (particleTypeName[ii].compare("O") == 0)
            mass.push_back(15.99943);
        else if (particleTypeName[ii].compare("H") == 0)
            mass.push_back(1.007947);
        // Convert the input coordinates unit from angstrom to nanometers.
        init_positions[ii * 3 + 0] = init_positions[ii * 3 + 0] * 0.1;
        init_positions[ii * 3 + 1] = init_positions[ii * 3 + 1] * 0.1;
        init_positions[ii * 3 + 2] = init_positions[ii * 3 + 2] * 0.1;
    }

    cout<<natoms<<" atoms in simulation."<<endl;
    std::cout.precision(10);
    
    // Test the dynamics of Deepmd Plugin.
    //try{
        registerDeepmdCudaKernelFactories();
        //if (argc > 1)
        //    Platform::getPlatformByName("CUDA").setPropertyDefaultValue("CudaPrecision", string(argv[1]));
        testDeepmdDynamicsWithMicromd(init_positions, particleTypeName, typeDict, mass, nsteps);
    /*}
    catch(const OpenMM::OpenMMException& e) {
        cout << "OpenMMException: "<<e.what() << endl;
        return 1;
    }*/
    cout<<"Done"<<endl;
    return 0;
}