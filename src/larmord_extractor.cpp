/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Aaron T. Frank
     
*/

#include "Molecule.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"
#include "LARMORD.hpp"
#include "Trajectory.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <math.h>
#include <algorithm>


void usage(){
  std::cerr << "====================================================" << std::endl;
  std::cerr << "LARMORD Chemical Shift Predictor v1.00" << std::endl;
  std::cerr << "(c) 2014 Aaron T. Frank and University of Michigan." << std::endl;
  std::cerr << "====================================================" << std::endl;
  std::cerr << "Usage:   larmord [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-csfile CSfile]" << std::endl;
  std::cerr << "         [-parmfile PARAMfile]" << std::endl;
  std::cerr << "         [-reffile REFfile]" << std::endl;
  std::cerr << "         [-accfile ACCfile]" << std::endl;
  std::cerr << "         [-corrfile CORRfile]" << std::endl;  
  std::cerr << "         [-cutoff CUToff]" << std::endl;
  std::cerr << "         [-printError]" << std::endl;
  std::cerr << "         [-residueSelection]" << std::endl;
  std::cerr << "         [-nucleusSelection]" << std::endl;
  std::cerr << "         [-residueBasedWeights]" << std::endl;  
  std::cerr << "         [-mismatchCheck]" << std::endl;  
  std::cerr << "         [-trj TRAJfile]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;  
  std::cerr << "         [-identification ID]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){
  int i;
  unsigned int f;
  unsigned int ainx;
  std::stringstream resid;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string fchemshift;
  std::string nucleus;
  std::string resname;
  std::string atomname;
  std::string key, residID;
  std::string identification;
  std::string residue_select;
  std::string nucleus_select;
  std::vector<int> selected_residues;
  std::vector<std::string> selected_nuclei;
  bool isResidue;
  bool isNucleus;
  bool mismatchCheck;
  std::vector<std::string> atomTypes;
  std::map<std::string, double> histo;
  
  std::vector<std::string> trajs;
  int start;
  int stop=std::numeric_limits<int>::max();
  int skip;
  bool startFlag=false;
  unsigned int itrj;
  std::ifstream trjin;
  Trajectory *ftrjin;
  unsigned int nframe;
  unsigned int process;
  unsigned int k;

  start=0;
  skip=0;
  nframe=0;
  process=0;
  
  int beta;
  double dist;
  double randcs;
  double expcs;
  double cutoff;
  

  std::vector<std::vector<double> > neighborDistances;

  Molecule *neighbormol;
  neighbormol=NULL;
  
  Atom *ai, *aj;  
  ai=NULL;
  aj=NULL;
  fchemshift="";
  identification="None";
  cutoff=99999.9;
  beta=-3;
  residue_select = "";
  nucleus_select = "";
  
 
  LARMORD *larm;
  larm=NULL;

  pdbs.clear();
  selected_residues.clear();
  selected_nuclei.clear();
  isResidue = true;
  mismatchCheck = false;
  atomTypes.clear();
  histo.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0)
    {
      usage();
    }
    else if (currArg.compare("-csfile") == 0 || currArg.compare("-c") == 0 )
    {
      currArg=argv[++i];
      fchemshift=currArg;
    }  
    else if (currArg.compare("-mismatchCheck") == 0 )
    {
				mismatchCheck=true;
    }    
    else if (currArg.compare("-residueSelection") == 0 )
    {
      currArg=argv[++i];    
			residue_select=currArg;			
    }      
    else if (currArg.compare("-nucleusSelection") == 0 )
    {
      currArg=argv[++i];    
			nucleus_select=currArg;
    }                
    else if (currArg.compare("-identification") == 0 || currArg.compare("-i") == 0 )
    {
      currArg=argv[++i];
      identification=currArg;
    }      
    else if (currArg.compare("-cutoff") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> cutoff;
    }
    else if (currArg.compare("-beta") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> beta;
    }
    else if (currArg.compare("-trj") == 0 || currArg.compare("-traj") == 0)
    {
      currArg=argv[++i];
      trajs.push_back(currArg);
    }
    else if (currArg.compare("-skip") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
      startFlag=true;
    }
    else if (currArg.compare("-stop") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }    
    else if (currArg.compare(0,1,"-") == 0)
    {
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }
  if (pdbs.size() == 0)
  {
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }
  
  Molecule *mol=NULL;

  if (residue_select.length() != 0 || nucleus_select.length() != 0){
		Molecule *cmol;
		cmol = Molecule::readPDB(pdbs.at(0));
		
		if (residue_select.length() != 0){
			cmol->select(residue_select);
			neighbormol=cmol->copy(true);	
			Residue *res;
			selected_residues.clear();
	
			for (unsigned int j=0; j< neighbormol->getResVecSize(); j++)
			{
				res=neighbormol->getResidue(j);
				selected_residues.push_back(res->getResId());
			}
		}

		if (nucleus_select.length() != 0){
			cmol->select(nucleus_select);
			neighbormol=cmol->copy(true);	
			Atom *atm;
			selected_nuclei.clear();			
			for (unsigned int j=0; j< neighbormol->getAtmVecSize(); j++)
			{
				atm=neighbormol->getAtom(j);
				if (std::find(selected_nuclei.begin(),selected_nuclei.end(),atm->getAtmName()) == selected_nuclei.end()){
					selected_nuclei.push_back(atm->getAtmName());
				}					
			}
		}
	}

  //There might be a cleaner way of doing this
  //but adding new atom types is easier/more obvious here.
  //Note that they are in alphanumeric order and sorted later!  
	atomTypes.push_back("CYS:C");
	atomTypes.push_back("CYS:CB");
	atomTypes.push_back("CYS:CA");
	atomTypes.push_back("CYS:O");
	atomTypes.push_back("CYS:N");
	atomTypes.push_back("CYS:SG");
	atomTypes.push_back("ASP:C");
	atomTypes.push_back("ASP:CB");
	atomTypes.push_back("ASP:CA");
	atomTypes.push_back("ASP:CG");
	atomTypes.push_back("ASP:O");
	atomTypes.push_back("ASP:N");
	atomTypes.push_back("ASP:OD1");
	atomTypes.push_back("ASP:OD2");
	atomTypes.push_back("SER:C");
	atomTypes.push_back("SER:OG");
	atomTypes.push_back("SER:CB");
	atomTypes.push_back("SER:CA");
	atomTypes.push_back("SER:O");
	atomTypes.push_back("SER:N");
	atomTypes.push_back("GLN:C");
	atomTypes.push_back("GLN:CB");
	atomTypes.push_back("GLN:CA");
	atomTypes.push_back("GLN:CG");
	atomTypes.push_back("GLN:O");
	atomTypes.push_back("GLN:N");
	atomTypes.push_back("GLN:CD");
	atomTypes.push_back("GLN:NE2");
	atomTypes.push_back("GLN:OE1");
	atomTypes.push_back("LYS:C");
	atomTypes.push_back("LYS:CB");
	atomTypes.push_back("LYS:CA");
	atomTypes.push_back("LYS:CG");
	atomTypes.push_back("LYS:CE");
	atomTypes.push_back("LYS:CD");
	atomTypes.push_back("LYS:NZ");
	atomTypes.push_back("LYS:O");
	atomTypes.push_back("LYS:N");
	atomTypes.push_back("ILE:C");
	atomTypes.push_back("ILE:CG2");
	atomTypes.push_back("ILE:CB");
	atomTypes.push_back("ILE:CA");
	atomTypes.push_back("ILE:O");
	atomTypes.push_back("ILE:N");
	atomTypes.push_back("ILE:CD1");
	atomTypes.push_back("ILE:CG1");
	atomTypes.push_back("ILE:CD");
	atomTypes.push_back("PRO:C");
	atomTypes.push_back("PRO:CB");
	atomTypes.push_back("PRO:CA");
	atomTypes.push_back("PRO:CG");
	atomTypes.push_back("PRO:O");
	atomTypes.push_back("PRO:N");
	atomTypes.push_back("PRO:CD");
	atomTypes.push_back("THR:C");
	atomTypes.push_back("THR:CB");
	atomTypes.push_back("THR:CA");
	atomTypes.push_back("THR:OG1");
	atomTypes.push_back("THR:O");
	atomTypes.push_back("THR:N");
	atomTypes.push_back("THR:CG2");
	atomTypes.push_back("PHE:C");
	atomTypes.push_back("PHE:CE1");
	atomTypes.push_back("PHE:CB");
	atomTypes.push_back("PHE:CA");
	atomTypes.push_back("PHE:CG");
	atomTypes.push_back("PHE:O");
	atomTypes.push_back("PHE:N");
	atomTypes.push_back("PHE:CZ");
	atomTypes.push_back("PHE:CD1");
	atomTypes.push_back("PHE:CD2");
	atomTypes.push_back("PHE:CE2");
	atomTypes.push_back("ALA:CB");
	atomTypes.push_back("ALA:C");
	atomTypes.push_back("ALA:CA");
	atomTypes.push_back("ALA:O");
	atomTypes.push_back("ALA:N");
	atomTypes.push_back("GLY:C");
	atomTypes.push_back("GLY:CA");
	atomTypes.push_back("GLY:O");
	atomTypes.push_back("GLY:N");
	atomTypes.push_back("HIS:C");
	atomTypes.push_back("HIS:CE1");
	atomTypes.push_back("HIS:CB");
	atomTypes.push_back("HIS:CA");
	atomTypes.push_back("HIS:CG");
	atomTypes.push_back("HIS:O");
	atomTypes.push_back("HIS:N");
	atomTypes.push_back("HIS:CD2");
	atomTypes.push_back("HIS:ND1");
	atomTypes.push_back("HIS:NE2");
	atomTypes.push_back("GLU:C");
	atomTypes.push_back("GLU:CB");
	atomTypes.push_back("GLU:CA");
	atomTypes.push_back("GLU:CG");
	atomTypes.push_back("GLU:O");
	atomTypes.push_back("GLU:N");
	atomTypes.push_back("GLU:OE2");
	atomTypes.push_back("GLU:CD");
	atomTypes.push_back("GLU:OE1");
	atomTypes.push_back("LEU:C");
	atomTypes.push_back("LEU:CB");
	atomTypes.push_back("LEU:CA");
	atomTypes.push_back("LEU:CG");
	atomTypes.push_back("LEU:O");
	atomTypes.push_back("LEU:N");
	atomTypes.push_back("LEU:CD1");
	atomTypes.push_back("LEU:CD2");
	atomTypes.push_back("ARG:C");
	atomTypes.push_back("ARG:CB");
	atomTypes.push_back("ARG:CA");
	atomTypes.push_back("ARG:CG");
	atomTypes.push_back("ARG:NE");
	atomTypes.push_back("ARG:O");
	atomTypes.push_back("ARG:CD");
	atomTypes.push_back("ARG:CZ");
	atomTypes.push_back("ARG:NH1");
	atomTypes.push_back("ARG:NH2");
	atomTypes.push_back("ARG:N");
	atomTypes.push_back("TRP:C");
	atomTypes.push_back("TRP:CD1");
	atomTypes.push_back("TRP:CZ2");
	atomTypes.push_back("TRP:CB");
	atomTypes.push_back("TRP:CA");
	atomTypes.push_back("TRP:CG");
	atomTypes.push_back("TRP:O");
	atomTypes.push_back("TRP:N");
	atomTypes.push_back("TRP:CH2");
	atomTypes.push_back("TRP:CE3");
	atomTypes.push_back("TRP:CE2");
	atomTypes.push_back("TRP:CD2");
	atomTypes.push_back("TRP:CZ3");
	atomTypes.push_back("TRP:NE1");
	atomTypes.push_back("VAL:C");
	atomTypes.push_back("VAL:CB");
	atomTypes.push_back("VAL:CA");
	atomTypes.push_back("VAL:O");
	atomTypes.push_back("VAL:N");
	atomTypes.push_back("VAL:CG1");
	atomTypes.push_back("VAL:CG2");
	atomTypes.push_back("ASN:C");
	atomTypes.push_back("ASN:CB");
	atomTypes.push_back("ASN:CA");
	atomTypes.push_back("ASN:CG");
	atomTypes.push_back("ASN:O");
	atomTypes.push_back("ASN:N");
	atomTypes.push_back("ASN:OD1");
	atomTypes.push_back("ASN:OT1");
	atomTypes.push_back("ASN:ND2");
	atomTypes.push_back("TYR:C");
	atomTypes.push_back("TYR:CE1");
	atomTypes.push_back("TYR:OH");
	atomTypes.push_back("TYR:CB");
	atomTypes.push_back("TYR:CA");
	atomTypes.push_back("TYR:CG");
	atomTypes.push_back("TYR:O");
	atomTypes.push_back("TYR:N");
	atomTypes.push_back("TYR:CZ");
	atomTypes.push_back("TYR:CD1");
	atomTypes.push_back("TYR:CD2");
	atomTypes.push_back("TYR:CE2");
	atomTypes.push_back("MET:C");
	atomTypes.push_back("MET:CB");
	atomTypes.push_back("MET:CA");
	atomTypes.push_back("MET:CG");
	atomTypes.push_back("MET:CE");
	atomTypes.push_back("MET:N");
	atomTypes.push_back("MET:O");
	atomTypes.push_back("MET:SD");
	
  
  if (trajs.size() > 0)
  {
    if (pdbs.size() > 1)
    {
      std::cerr << std::endl << "Warning: Only the first PDB structure is used for trajectory analysis" << std::endl << std::endl;
    }
    /* Trajectory analysis */
    /* instantiate LARMORD */
    mol=Molecule::readPDB(pdbs.at(0));
    mol->selAll();
    larm = new LARMORD(mol,fchemshift,"","","","",false,false,false,true);
        
    /* Process trajectories */
    for (itrj=0; itrj< trajs.size(); itrj++)
    {
      trjin.open(trajs.at(itrj).c_str(), std::ios::binary);
      if (trjin.is_open())
      {
        ftrjin=new Trajectory;
        ftrjin->setMolecule(mol);
        if (ftrjin->findFormat(trjin) == true)
        {
          ftrjin->readHeader(trjin);
          if (skip > 0 && startFlag == false)
          {
            start=skip;
          }
          /* Loop through desired frames */
          for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip)
          {
            if( ftrjin->readFrame(trjin, i) == false)
            {
              std::cerr << "Warning: EOF found before the next frame could be read" << std::endl;
              break;
            }
            nframe++;
            /* get distance matrix */
            mol->assignAtmInx();
            mol->selAll();
            Analyze::pairwiseDistance(mol, neighborDistances);
            
            mol->select(":.HEAVY");
            neighbormol= mol->copy();
         
						for (unsigned int j=0; j< mol->getAtmVecSize(); j++)
						{
							ai = mol->getAtom(j);
							nucleus = ai->getAtmName();
							if (larm->getShiftAtom(nucleus)==true)
							{
								resname = ai->getResName();
								resid << ai->getResId();
								residID = resid.str();
								key = resid.str()+":"+resname+":"+nucleus;
								resid.str("");            
								expcs = larm->getExperimentalCS(key);
								//std::cout << key << expcs << std::endl;
								isResidue = std::find(selected_residues.begin(),selected_residues.end(),ai->getResId())!= selected_residues.end();
								isNucleus = std::find(selected_nuclei.begin(),selected_nuclei.end(),ai->getAtmName()) != selected_nuclei.end();
								if( (fchemshift.length() == 0 || expcs != 0.0) && (selected_residues.size() == 0 || isResidue) && (selected_nuclei.size() == 0 || isNucleus) )
								{
									key = resname+":"+nucleus;
									randcs = larm->getRandomShift(key);
									histo.clear();        																	
									if(true)
									{
										ainx = ai->getAtmInx();
										for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++)
										{
											aj = neighbormol->getAtom(l);
											dist = neighborDistances.at(ainx).at(aj->getAtmInx());
								
											if(ai!=aj && dist < cutoff)
											{
												key = aj->getResName()+":"+aj->getAtmName();												
												if (histo.find(key) != histo.end()){
													histo.at(key)+=pow(dist,beta);
												}
												else{
													histo.insert(std::pair<std::string, double>(key, pow(dist,beta)));
												}
											}
										}
									}
									process++;
									// print header
									if(process==1)
									{
										std:: cout << "ID resname resid nucleus expCS ";
										for (k=0; k< atomTypes.size(); k++){
											key=atomTypes.at(k);
											std::cout << key;
											if (k != atomTypes.size()-1 ){
												std::cout << " ";
											}
										}
										std::cout << std::endl;  
									} 
									// print histogram
									std:: cout << identification << " " << resname << " " << residID << " " <<  nucleus << " " << expcs << " "; 
									for (k=0; k< atomTypes.size(); k++)
									{
										key=atomTypes.at(k);
										if (histo.find(key) != histo.end()){
											std::cout << histo.at(key);
										}
										else{
											std::cout << "0";
										}
										if (j != atomTypes.size()-1){
											std::cout << " ";
										}
									}
									std::cout << std::endl;															
								}
							}        
						}
            
          }
        }
        else
        {
          std::cerr << "Warning: Skipping unknown trajectory format \"";
          std::cerr << trajs.at(itrj) << "\"" << std::endl;
        }
        if (ftrjin != NULL){
          delete ftrjin;
        }
      }
      trjin.close();
    }
  }
  else 
  { 
    /* instantiate LARMORD */
    for (f=0; f< pdbs.size(); f++)
    {  
      mol=Molecule::readPDB(pdbs.at(f));
      larm = new LARMORD(mol,fchemshift,"","","","",false,false,false,true);
      
      //std::cerr << "Processing file \"" << pdbs.at(f) << "..." << std::endl;
      /* get distance matrix */
      mol->assignAtmInx();
      mol->selAll();
      Analyze::pairwiseDistance(mol, neighborDistances);

      mol->select(":.HEAVY");
      neighbormol= mol->copy();
            
      for (unsigned int j=0; j< mol->getAtmVecSize(); j++)
      {
        ai = mol->getAtom(j);
        nucleus = ai->getAtmName();
        if (larm->getShiftAtom(nucleus)==true)
        {
          resname = ai->getResName();
          resid << ai->getResId();
          residID = resid.str();
          key = resid.str()+":"+resname+":"+nucleus;
          resid.str("");            
          expcs = larm->getExperimentalCS(key);
          //std::cout << key << expcs << std::endl;
          isResidue = std::find(selected_residues.begin(),selected_residues.end(),ai->getResId())!= selected_residues.end();
          isNucleus = std::find(selected_nuclei.begin(),selected_nuclei.end(),ai->getAtmName()) != selected_nuclei.end();
          if( (fchemshift.length() == 0 || expcs != 0.0) && (selected_residues.size() == 0 || isResidue) && (selected_nuclei.size() == 0 || isNucleus) )
          {
						key = resname+":"+nucleus;
						randcs = larm->getRandomShift(key);
						histo.clear();        																	
						if(true)
						{
							ainx = ai->getAtmInx();
							for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++)
							{
								aj = neighbormol->getAtom(l);
								dist = neighborDistances.at(ainx).at(aj->getAtmInx());
								
								if(ai!=aj && dist < cutoff)
								{
									key = aj->getResName()+":"+aj->getAtmName();												
									if (histo.find(key) != histo.end()){
										histo.at(key)+=pow(dist,beta);
									}
									else{
										histo.insert(std::pair<std::string, double>(key, pow(dist,beta)));
									}
								}
							}
						}
						process++;
						// print header
						if(process==1)
						{
							std:: cout << "ID resname resid nucleus expCS ";
							for (k=0; k< atomTypes.size(); k++){
								key=atomTypes.at(k);
								std::cout << key;
								if (k != atomTypes.size()-1 ){
									std::cout << " ";
								}
							}
							std::cout << std::endl;  
						} 
						// print histogram
						std:: cout << identification << " " << resname << " " << residID << " " <<  nucleus << " " << expcs << " "; 
						for (k=0; k< atomTypes.size(); k++)
						{
							key=atomTypes.at(k);
							if (histo.find(key) != histo.end()){
								std::cout << histo.at(key);
							}
							else{
								std::cout << "0";
							}
							if (j != atomTypes.size()-1){
								std::cout << " ";
							}
						}
						std::cout << std::endl;															
          }
        }        
      }
      delete mol;
    }
  }  
  if(larm!=NULL)
    delete larm;
  return 0;
}
