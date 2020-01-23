/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
    ooo                        | 
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef LEMONADE_ANALYZER_ANALYZERBONDLENGTH
#define LEMONADE_ANALYZER_ANALYZERBONDLENGTH

#include <string>
#include <iostream>

#include <LeMonADE/Version.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/DistanceCalculation.h>


/***********************************************************************/
/**
 * @file
 *
 * @class AnalyzerBondLength
 * 
 * @brief Analyzer calculating Radius of Gyration of dendrimer and averaged Radius of Gyration of all solvent chains
 * 
 * @details 
 * 
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 **/
template <class IngredientsType> 
class AnalyzerBondLength: public AbstractAnalyzer
{
public:
  //! constructor call transferring ingredients and the result file name
  AnalyzerBondLength(const IngredientsType& ing, const std::string& filename_);

  //!Standard destructor closing the file stream.
  virtual ~AnalyzerBondLength();

  //! initialise the groups to calculate densities and define the middle points
  virtual void initialize();
  
  //! do the calculation for the actual mcs
  virtual bool execute();
  
  //! performs all the result to file operations
  virtual void cleanup();
  
  //! Returns the result filename used in this class
  std::string getFilename(){return filename;}
  
  //! getter for initialised bool
  bool getIsInitialized(){ return isInitialized;}
  
  //! getter for numiber of executions
  int32_t getNumExec(){ return numExec;}
  
  //! get size of result vector
  size_t getSizeResultVector(){ return results_bond_length.size();}
  
private:
  //! Storage for data that are processed to file (mostly Ingredients).
  const IngredientsType& ingredients;
  
  //! The filename to store the results in
  std::string filename;
  
  //! bool to check if groups are initilized
  bool isInitialized;
  
  //! number of executions
  int32_t numExec;
  
  //! container for results
  std::vector< std::vector < double > > results_bond_length;
  
  std::vector<double>averaged_bond_length;
  
};

/***********************************************************************/
//constructor
/***********************************************************************/
/**
 * @details It initialize all internal values and passes all information to the corresponding classes.
 *
 * @param filename Name of the file to write-out
 * @param ingredients_ Class holding all information of the system (mainly Ingredients )
 *
 * @todo Did I understand that correctly that the first !mcs is read not the last in the file?
 */
template <class IngredientsType>
AnalyzerBondLength<IngredientsType>::AnalyzerBondLength(const IngredientsType& ingredients_, const std::string& filename_)
    :ingredients(ingredients_),
    filename(filename_),isInitialized(false), numExec(0), 
    results_bond_length(3,std::vector<double>())
{}

/***********************************************************************/
//destructor
/***********************************************************************/
/**
 * @details it frees the memory from the write objects, which were registered
 * as pointers with registerWrite
 */
template<class IngredientsType>
AnalyzerBondLength<IngredientsType>::~AnalyzerBondLength()
{}


/***********************************************************************/
//void initialize
/***********************************************************************/
/**
 * @details find the dendrimer in the box and the linear chains
 *
 */
template <class IngredientsType>
void AnalyzerBondLength<IngredientsType>::initialize(){
  std::cout <<"initialize AnalyzerBondLength"<<std::endl;
  
  isInitialized=true;
  
  execute();
}

/***********************************************************************/
//bool execute
/***********************************************************************/
/**
 * @details calculate overla with the densities of the current mcs
 *
 * @return True if everthing is alrigth. False if something goes wrong.
 */
template <class IngredientsType>
bool AnalyzerBondLength<IngredientsType>::execute(){
  // check if initialise was called before
  if(isInitialized==false) initialize();
  
  // count number of execution
  numExec++;
  
  // throw away first 10 configs for equlibration
  //if(numExec>10){
    //std::cout << "execute frame number "<<numExec<<std::endl;
    double actual_mean_squared_bondlength(0.0);
    int32_t number_of_counted_bonds(0);
    
    // collect all bonds (twice)
    for(uint32_t i=0; i<ingredients.getMolecules().size();i++){
      for(uint32_t j=0;j<ingredients.getMolecules().getNumLinks(i);j++){
	      double length( (ingredients.getMolecules()[i]-ingredients.getMolecules()[ingredients.getMolecules().getNeighborIdx(i,j)]).getLength() );
	      actual_mean_squared_bondlength+=length*length;
	
	      number_of_counted_bonds++;
      }
    }
    
    //check if there were some bonds
    if(number_of_counted_bonds==0)
      throw std::runtime_error("no bonds to analyze!");
    
    // add actual bond to averaging container
    averaged_bond_length.push_back(actual_mean_squared_bondlength/((double)number_of_counted_bonds));
    
    //fill the results vector: 
    // mcs  actual bondlength   averaged bondlength
    results_bond_length.at(0).push_back(ingredients.getMolecules().getAge()); 
    results_bond_length.at(1).push_back(actual_mean_squared_bondlength/((double)number_of_counted_bonds));
    double sum(0.0);
    for(size_t i=0;i<averaged_bond_length.size();i++){
      sum+=averaged_bond_length.at(i);
    }
    sum/=(averaged_bond_length.size());
    results_bond_length.at(2).push_back(sum);
  //}
  
  return true;
}

/***********************************************************************/
//void cleanup
/***********************************************************************/
/**
 * @details write out results of Overlap Analyzer to file
 */
template <class IngredientsType>
void AnalyzerBondLength<IngredientsType>::cleanup(){
  //write out files
  std::cout << "clean up AnalyzerBondLength"<<std::endl;
  // use the result formatting tools to write out the results
  ResultFormattingTools::writeResultFile(filename, ingredients, results_bond_length,"mean squared bond length\nmcs\tBL frame \tBL averaged");
}

#endif //LEMONADE_ANALYZER_ANALYZERBONDLENGTH
