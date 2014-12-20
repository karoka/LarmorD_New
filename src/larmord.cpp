/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#include "Molecule.hpp"
#include "Atom.hpp"
#include "Analyze.hpp"
#include "LARMORD.hpp"
#include "Trajectory.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <math.h>


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
  std::cerr << "         [-cutoff CUToff]" << std::endl;
  std::cerr << "         [-printError]" << std::endl;
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
  unsigned int counter;
  std::stringstream resid;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string fchemshift;
  std::string fparmfile;
  std::string freffile;
  std::string faccfile;
  std::string nucleus;
  std::string resname;
  std::string atomname;
  std::string key;
  std::string identification;
  bool print_error;

  std::vector<std::string> trajs;
  int start;
  int stop=std::numeric_limits<int>::max();
  int skip;
  bool startFlag=false;
  unsigned int itrj;
  std::ifstream trjin;
  Trajectory *ftrjin;
  unsigned int nframe;

  start=0;
  skip=0;
  nframe=0;
  
  std::vector<double> alpha;
  std::vector<int> beta;
  double dist;
  double cspred;
  double randcs;
  double expcs;
  double cutoff;
  double error, error_mae, error_rmse, error_wmae, error_wrmse, weight;

  std::vector<std::vector<double> > neighborDistances;

  Molecule *neighbormol;
  neighbormol=NULL;
  
  Atom *ai, *aj;  
  ai=NULL;
  aj=NULL;
  fchemshift="";
  fparmfile="";
  freffile="";  
  faccfile="";  
  identification="None";
  cutoff=99999.9;
  print_error = false;
  
  LARMORD *larm;
  larm=NULL;

  pdbs.clear();

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
    else if (currArg.compare("-parmfile") == 0 || currArg.compare("-p") == 0 )
    {
      currArg=argv[++i];
      fparmfile=currArg;
    }
    else if (currArg.compare("-reffile") == 0 || currArg.compare("-r") == 0 )
    {
      currArg=argv[++i];
      freffile=currArg;
    }
    else if (currArg.compare("-printError") == 0 )
    {
			if(fchemshift.length() > 0)
			{
					print_error=true;
			} else 
			{
					std::cerr << "Warning: -printError ignored" << std::endl;
					std::cerr << "Warning: -printError is valid only when a chemical shifts data file is suppled via --csfile " << std::endl;
					print_error=false;
			}      
    }    
    else if (currArg.compare("-accfile") == 0 || currArg.compare("-a") == 0 )
    {
      currArg=argv[++i];
      faccfile=currArg;
    }
    else if (currArg.compare("-identification") == 0 || currArg.compare("-i") == 0 )
    {
      currArg=argv[++i];
      identification=currArg;
    }      
    else if (currArg.compare("-trj") == 0 || currArg.compare("-traj") == 0)
    {
      currArg=argv[++i];
      trajs.push_back(currArg);
    }
    else if (currArg.compare("-cutoff") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> cutoff;
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
    larm = new LARMORD(mol,fchemshift,fparmfile,freffile,faccfile);
    
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

						if (print_error)
						{
							error_mae = 0.0;
							error_rmse = 0.0;
							error_wmae = 0.0;
							error_wrmse = 0.0;
							counter = 0;
						}
            
            for (unsigned int j=0; j< mol->getAtmVecSize(); j++)
            {
              ai = mol->getAtom(j);
              nucleus = ai->getAtmName();
              if (larm->getShiftAtom(nucleus)==true)
              {
                resname = ai->getResName();
                resid << ai->getResId();
                key = resid.str()+":"+nucleus;
                resid.str("");
                cspred = 0.0;
                expcs = larm->getExperimentalCS(key);
                //std::cerr << " I am here " << key << " " << expcs << std::endl;
                if(fchemshift.length() == 0 || expcs != 0.0)
                {
									key = resname+":"+nucleus;
									randcs = larm->getRandomShift(key);
									if(randcs > 0)
									{
										ainx = ai->getAtmInx();
										for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++)
										{
											aj = neighbormol->getAtom(l);
											if(ai!=aj){
												resname = aj->getResName();
												atomname = aj->getAtmName();
												alpha = larm->getAlpha(nucleus+":"+aj->getResName()+":"+aj->getAtmName());
												beta = larm->getBeta(nucleus+":"+aj->getResName()+":"+aj->getAtmName());
												dist = neighborDistances.at(ainx).at(aj->getAtmInx());
												if (dist < cutoff){
													for (unsigned int m = 0; m < alpha.size(); m++)
													{
														cspred = cspred + alpha.at(m)*pow(dist,beta.at(m));
													}
												}                          
											}
										}
										cspred = cspred + randcs;
										if(print_error)
										{
											weight = larm->getAccuracyWeight(nucleus);
											cspred = cspred + randcs;
											error = cspred - expcs;
											error_mae += fabs(error);
											error_rmse += (error*error);
											error_wmae += weight*fabs(error);
											error_wrmse += weight*weight*(error*error);
											counter++; 
										} 
										else 
										{                               
											std::cout << 0 << " " << nframe << " " << ai->getResId() << " " << ai->getResName() << " " << nucleus << " " << cspred << " " << expcs << " " <<  randcs << " " << identification << std::endl;
										}                      
									}
                }
              }
            }
						if (print_error)
						{
							std::cout << 0 << " " << nframe << " " << error_mae/counter << " " << sqrt(error_rmse/counter) << " " << error_wmae/counter << " " <<  sqrt(error_wrmse/counter) << " " << identification << std::endl;
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
      larm = new LARMORD(mol,fchemshift,fparmfile,freffile,faccfile);
      //std::cerr << "Processing file \"" << pdbs.at(f) << "..." << std::endl;
      /* get distance matrix */
      mol->assignAtmInx();
      mol->selAll();
      Analyze::pairwiseDistance(mol, neighborDistances);

      mol->select(":.HEAVY");
      neighbormol= mol->copy();
      
      if (print_error)
      {
				error_mae = 0.0;
				error_rmse = 0.0;
				error_wmae = 0.0;
				error_wrmse = 0.0;
				counter = 0;
      }
      
      for (unsigned int j=0; j< mol->getAtmVecSize(); j++)
      {
        ai = mol->getAtom(j);
        nucleus = ai->getAtmName();
        if (larm->getShiftAtom(nucleus)==true)
        {
          resname = ai->getResName();
          resid << ai->getResId();
          key = resid.str()+":"+nucleus;
          resid.str("");
          cspred = 0.0;
          expcs = larm->getExperimentalCS(key);
          //std::cerr << " I am here " << key << " " << expcs << std::endl;
          if(fchemshift.length() == 0 || expcs != 0.0)
          {
						key = resname+":"+nucleus;
						randcs = larm->getRandomShift(key);
						if(randcs > 0)
						{
							ainx = ai->getAtmInx();
							for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++){
								aj = neighbormol->getAtom(l);
								if(ai!=aj){
									resname = aj->getResName();
									atomname = aj->getAtmName();
									alpha = larm->getAlpha(nucleus+":"+aj->getResName()+":"+aj->getAtmName());
									beta = larm->getBeta(nucleus+":"+aj->getResName()+":"+aj->getAtmName());                    
									dist = neighborDistances.at(ainx).at(aj->getAtmInx());
									if (dist < cutoff)
									{
										for (unsigned int m = 0; m < alpha.size(); m++)
										{
											//std::cout << nucleus+":"+aj->getResName()+":"+aj->getAtmName() << " " << alpha.at(m) << " " << beta.at(m) << " " << randcs << std::endl;
											cspred = cspred + alpha.at(m)*pow(dist,beta.at(m));
										}
									}
								}
							}							
							if(print_error)
							{
								weight = larm->getAccuracyWeight(nucleus);
								cspred = cspred + randcs;
								error = cspred - expcs;
								error_mae += fabs(error);
								error_rmse += (error*error);
								error_wmae += weight*fabs(error);
								error_wrmse += weight*weight*(error*error);
								counter++; 
							} 
							else 
							{                               
								std::cout << 0 << " " << f+1 << " " << ai->getResId() << " " << ai->getResName() << " " << nucleus << " " << cspred << " " << expcs << " " <<  randcs << " " << identification << std::endl;
							}
						}
          }
        }
      }
      if (print_error)
      {
      	std::cout << 0 << " " << f+1 << " " << error_mae/counter << " " << sqrt(error_rmse/counter) << " " << error_wmae/counter << " " <<  sqrt(error_wrmse/counter) << " " << identification << std::endl;
      }
      delete mol;
    }
  }  
  if(larm!=NULL)
    delete larm;
  return 0;
}
