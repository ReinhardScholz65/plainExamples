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

#ifndef LEMONADE_ANALYZER_PAIR_DISTRIBUTION
#define LEMONADE_ANALYZER_PAIR_DISTRIBUTION

#include<vector>

#include <LeMonADE/analyzer/AbstractAnalyzer.h>

/***********************************************************************/
/**
 * @file
 *
 * @class AnalyzerPairDistribution
 **/
template <class IngredientsType> 
class AnalyzerPairDistribution: public AbstractAnalyzer
{
public:
  //AnalyzerPairDistribution(const IngredientsType& ing, int chainlength);
  AnalyzerPairDistribution(const IngredientsType& ing);

  virtual ~AnalyzerPairDistribution();

  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();
  
private:
  //! const reference to ingredients
  const IngredientsType& ingredients;

  //! bool to check if analyzer was initialized
  bool wasInitialized;
  
};

/**
 * @details constructur passing ingredients and your parameters to the system
 *
 * @param myParameter_ a parameter to be passed to the analyzer 
 * @param ingredients_ Class holding all information of the system (mainly Ingredients )
 */
template <class IngredientsType>
AnalyzerPairDistribution<IngredientsType>::AnalyzerPairDistribution(const IngredientsType& ingredients_):ingredients(ingredients_), wasInitialized(false){}

/**
 * @details default (empty) destructor
 * 
 */
template<class IngredientsType>
AnalyzerPairDistribution<IngredientsType>::~AnalyzerPairDistribution()
{}

/**
 * @details The initialize function handles the new systems information before execution.
 *
 */
template <class IngredientsType>
void AnalyzerPairDistribution<IngredientsType>::initialize(){
    std::cout << "initialize AnalyzerPairDistribution" << std::endl;

    // std::cout << "you passed the parameter "<< myParameter <<std::endl;

    wasInitialized=true;
    
    // by default, initialise calls execute
    execute();
}

/**
 * @details Execution of your task.
 * 
 */
template <class IngredientsType>
bool AnalyzerPairDistribution<IngredientsType>::execute(){
  if(wasInitialized==false) initialize();
  std::cout <<"execute AnalyzerPairDistribution"<<std::endl;

  int numberMonomers;
  numberMonomers = ingredients.getMolecules().size();
  std::cout << "number of monomers: " << numberMonomers << std::endl;
 
  int i2; 
  int j; 
  int k; 
  int ix; 
  int iy; 
  int iz;

  VectorInt3 position_j; 
  VectorInt3 position_k;
  VectorInt3 delta_jk;
  
  std::vector<int> g2(1025);  // or directly: g2(1025,0);
  // VectorInt g2(1025);

  for(i2 = 0; i2 < 1024; ++i2) {
    g2[i2] = 0;
    // std::cout << "pair distribution(" << i2 << ") = " << g2[i2] << std::endl;
  }
  //
  for(j = 0; j < numberMonomers; j++) {
    position_j = ingredients.getMolecules()[j].getVector3D();
    for(k = 0; k < numberMonomers; k++) {
      position_k = ingredients.getMolecules()[k].getVector3D();
      delta_jk = position_j - position_k;
      ix = delta_jk[0]; 
      iy = delta_jk[1];       
      iz = delta_jk[2];
      i2 = ix * ix + iy * iy +iz *iz;
      //std::cout << "pair(" << j << "," << k << ") : (" << position_j << ") and (" << position_k << "), delta = (" 
      //          << delta_jk << "), delta**2 = " << i2 << std::endl;
      if(i2 <= 1024) {
          g2[i2] = g2[i2] + 1;
      }
    }
  }

for(i2 = 0; i2 < 200; ++i2) {
    if(g2[i2]!=0) {
    // std::cout << "pair distribution(" << i2 << ") = " << g2[i2] << std::endl;
    std::cout << i2 << " " << g2[i2] << std::endl;
    }
  }
 
  
  return true;
}

/**
 * @details write out results of Overlap Analyzer to file
 */
template <class IngredientsType>
void AnalyzerPairDistribution<IngredientsType>::cleanup(){
  //example of write out files:
  // you need to include: 
  //    #include <LeMonADE/utility/ResultFormattingTools.h>
  // and run
  // ResultFormattingTools::writeResultFile(filename, ingredients, some std::vector<std::vector<TYPE> >,"some information about the data in the file\nColumnsA\tColumnsB");
}

#endif //LEMONADE_ANALYZER_ANALYZERMYRGDENDRIMER
