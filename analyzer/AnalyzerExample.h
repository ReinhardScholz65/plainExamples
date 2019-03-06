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

#ifndef LEMONADE_ANALYZER_ANALYZERMYRGDENDRIMER
#define LEMONADE_ANALYZER_ANALYZERMYRGDENDRIMER

#include <LeMonADE/analyzer/AbstractAnalyzer.h>

/***********************************************************************/
/**
 * @file
 *
 * @class AnalyzerExample
 **/
template <class IngredientsType> 
class AnalyzerExample: public AbstractAnalyzer
{
public:
  AnalyzerExample(const IngredientsType& ing, const int32_t myParameter_);

  virtual ~AnalyzerExample();

  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();
  
private:
  //! const reference to ingredients
  const IngredientsType& ingredients;
  
  //! a paramter you might have to use
  const uint32_t myParameter;

  //! bool to check if analyzer was initialized
  bool wasInitialized;
  
};

/**
 * @details constructur passing ingredients and your paramters to the system
 *
 * @param myParameter_ a paramter to be passed to the analyzer 
 * @param ingredients_ Class holding all information of the system (mainly Ingredients )
 */
template <class IngredientsType>
AnalyzerExample<IngredientsType>::AnalyzerExample(const IngredientsType& ingredients_, const int32_t myParameter_):
    ingredients(ingredients_),myParameter(myParameter_), wasInitialized(false)
{}

/**
 * @details default (empty) destructor
 * 
 */
template<class IngredientsType>
AnalyzerExample<IngredientsType>::~AnalyzerExample()
{}

/**
 * @details The initialize function handles the new systems information before execution.
 *
 */
template <class IngredientsType>
void AnalyzerExample<IngredientsType>::initialize(){
    std::cout << "initialize AnalyzerExample" << std::endl;

    std::cout << "you passed the paramter "<< myParameter <<std::endl;

    wasInitialized=true;
    
    // by default, initialise calls execute
    execute();
}

/**
 * @details Execution of your task.
 * 
 */
template <class IngredientsType>
bool AnalyzerExample<IngredientsType>::execute(){
  if(wasInitialized==false) initialize();
  std::cout <<"execute AnalyzerExample"<<std::endl;
  
  return true;
}

/**
 * @details write out results of Overlap Analyzer to file
 */
template <class IngredientsType>
void AnalyzerExample<IngredientsType>::cleanup(){
  //example of write out files:
  // you need to include: 
  //    #include <LeMonADE/utility/ResultFormattingTools.h>
  // and run
  // ResultFormattingTools::writeResultFile(filename, ingredients, some std::vector<std::vector<TYPE> >,"some information about the data in the file\nColumnsA\tColumnsB");
}

#endif //LEMONADE_ANALYZER_ANALYZERMYRGDENDRIMER
