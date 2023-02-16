#include <iostream>
#include "openmm/OpenMMException.h"
#include "openmm/Context.h"
#include "openmm/reference/SimTKOpenMMRealType.h"


using namespace std;
using namespace OpenMM;


namespace DeepmdPlugin {

class OPENMM_EXPORT_DEEPMD Topology;

class OPENMM_EXPORT_DEEPMD Atom{
public:
    Atom(){
        name = "";
        element = "";
        index = -1;
        ResIndex = -1;
        id = -1;
        DisWithZinc = 0.0;
        position = Vec3();
    }
    // Initialize the instance with input only.
    Atom::Atom(const int ResIndex, string name, string element, int index, int id) : name(name), element(element), index(index), ResIndex(ResIndex), id(id){
        position = Vec3();
    }
    ~Atom(){};
    const string getElement() const{
        return element;
    }
    const int getIndex() const{
        return index;
    }
    const int getResIndex() const {
        return ResIndex;
    }

    void setPosition(Vec3 pos){
        position = pos;
    }

    std::string name;
    std::string element;
    int index;
    int ResIndex;
    int id; // atom serial number in .pdb file actually.

    Vec3 position;
};

class OPENMM_EXPORT_DEEPMD Residue{
public:
    Residue(){
        atoms = vector<Atom>();
        ResName = "";
        ResIndex = -1;
        ChainIndex = -1;
        ResId = -1;
    }

    Residue(const int chainIndex, std::string name, int ResIndex, int ResId): ChainIndex(chainIndex), ResName(name), ResIndex(ResIndex), ResId(ResId){
        atoms = vector<Atom>();
    }
    
    ~Residue(){};

    const vector<Atom> getAtoms() const {
        return atoms;
    }
    const vector<int> getAtomsIndex() const {
        return atoms_index;
    }
    const int getResIndex() const{
        return ResIndex;
    }
    const int getChainIndex() const{
        return ChainIndex;
    }
    std::vector<Atom> atoms;
    vector<int> atoms_index;
    std::string ResName;
private:
    int ResIndex;
    int ChainIndex;
    int ResId; // Residue serial number in .pdb file.
};

class OPENMM_EXPORT_DEEPMD Chain{
public:
    Chain(){
        residues = vector<Residue>();
        ChainIndex = -1;
        topology = NULL;
        ChainId = -1;
    }

    Chain(int ChainIndex, const Topology* top, int ChainId) : ChainIndex(ChainIndex), topology(top), ChainId(ChainId){
        residues = vector<Residue>();
    }
    
    ~Chain(){};

    const vector<Residue> getResidues() const {
        return residues;
    }
    const vector<int> getResiduesIndex() const {
        return residues_index;
    }
    const int getChainIndex() const {
        return ChainIndex;
    }
    std::vector<Residue> residues;
    vector<int> residues_index;
private:
    int ChainIndex;
    const Topology *topology;
    int ChainId;    
};

bool operator==(const Atom& at1, const Atom& at2){
    return at1.getIndex() == at2.getIndex();
}
bool operator==(const Residue& res1, const Residue& res2){
     return res1.getResIndex() == res2.getResIndex();
}
bool operator==(const Chain& chain1, const Chain& chain2){
    return chain1.getChainIndex() == chain2.getChainIndex();
}

class OPENMM_EXPORT_DEEPMD Topology{
public:
    Topology(){
        numParticles = 0;
        chains = map<int, Chain>();
        atoms = map<int, Atom>();
        residues = map<int, Residue>();
        bonds = vector<vector<int>>();
    }

    ~Topology(){};
    
    int addChain(int chainIndex, int chainId){
        int oldSize = chains.size();
        chains[chainIndex] = Chain(chainIndex, this, chainId);

        int newSize = chains.size();
        if((newSize - oldSize)!=1){
            throw OpenMMException("DeepmdPlugin::Topology: add Chain failed.");
        }
        // Return the index of this chain in chains. So that you can call this new added chain with 
        // top.chains[idx]. where idx is the return value of this function. 
        return newSize - 1;
    }

    int addResidue(int chainIndex, string ResName, int ResIndex, int ResId){
        int oldSize = residues.size();

        Residue res = Residue(chainIndex, ResName, ResIndex, ResId);
        // Put the new added residues into topology residues.
        residues[ResIndex] = res;
        chains[chainIndex].residues.push_back(res);
        chains[chainIndex].residues_index.push_back(ResIndex);
        int newSize = residues.size();
        if((newSize - oldSize)!=1){
            throw OpenMMException("DeepmdPlugin::Topology: add residues failed.");
        }
        return newSize - 1;
    }

    int addAtom(int resIndex, string AtomName, string AtomElement, int atomIndex, int atomId){
        int oldSize = atoms.size();
        Atom at = Atom(resIndex, AtomName, AtomElement, atomIndex, atomId);
        atoms[atomIndex] = at;
        residues[resIndex].atoms.push_back(at);
        residues[resIndex].atoms_index.push_back(atomIndex);

        numParticles +=  1;
        int newSize = atoms.size();
        if ((newSize - oldSize)!=1){
            throw OpenMMException("DeepmdPlugin::Topology: add residues failed.");
        }
        return newSize - 1;
    }

    void addBond(int atomIndex1, int atomIndex2){
        if(atomIndex1 > bonds.size() || atomIndex2 > bonds.size()){
            throw OpenMMException("DeepmdPlugin::Topology: add bond failed. Atom index overflow.");
        }
        bonds[atomIndex1].push_back(atomIndex2);
        bonds[atomIndex2].push_back(atomIndex1);
    }

    const map<int, Atom>& getAtoms() const {
        return atoms;
    }
    const map<int, Residue>& getResidues() const{
        return residues;
    }
    const map<int, Chain>& getChains() const {
        return chains;
    }
    const vector<vector<int> >& getBonds() const {
        return bonds;
    }

private:
    int numParticles = 0;
    map<int, Chain> chains;
    map<int, Atom> atoms;
    map<int, Residue> residues;
    vector<vector<int> > bonds;
};

map<string, vector<int>> SearchAtomsInRegion(vector<Vec3> pos, vector<int> center_atoms, double radius, Topology* top, vector<string> atom_names, vector<int> &addOrNot){
    


}

}