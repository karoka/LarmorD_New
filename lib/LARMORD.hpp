/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef LARMORD_H
#define LARMORD_H

#include <string>
#include <vector>
#include <map>


//#include "Molecule.hpp"

/* Forward Declaration  (only valid for pointers and references) */
 class Molecule;

class LARMORD {
    private:
         std::map<std::string,bool> shiftAtoms;
         std::map<std::string,double> randomShifts;
         std::map<std::string, std::vector < double>  > alphas;
         std::map<std::string,std::vector < int> > betas;
         std::map<std::string,double> experimentalCS; 
         std::map<std::string,double> accuracy_weight; 
    public:
        LARMORD (Molecule *mol=NULL, const std::string fchemshift="",const std::string fparmfile="",const std::string freffile="",const std::string faccfile="");
        void initializeShiftAtoms();
        void initializeRandomShifts();
        void initializeAlpha();
        void initializeAccuracyWeight();
        bool getShiftAtom(const std::string &key);
        double getRandomShift(const std::string &key);
        std::vector<double> getAlpha (const std::string &key);
        std::vector<int> getBeta (const std::string &key);
        double getExperimentalCS(const std::string &key);
        double getAccuracyWeight(const std::string &key);
        int getNShiftAtoms();
        void renameRes(Molecule *mol);
        void loadCSFile(const std::string fchemshift);
        void loadParmFile(const std::string fparmfile);
        void loadRefFile(const std::string freffile);
        void loadAccFile(const std::string faccfile);
};
#endif
