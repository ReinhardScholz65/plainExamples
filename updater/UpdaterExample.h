/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
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

#ifndef LEMONADE_UPDATER_YVESEXAMPLE_H
#define LEMONADE_UPDATER_YVESEXAMPLE_H

/**
 * @file
 * @class UpdaterExample
 *
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>

template<class IngredientsType>
class UpdaterExample: public AbstractUpdater
{
  
public:
  UpdaterExample(IngredientsType& ingredients_, uint32_t myParameter_);
  
  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();
  
private:
  //! reference to ingredients
  IngredientsType& ingredients;
  
  //! a paramter you might have to use
  uint32_t myParameter;
  
  //! bool to ensure updater runs only ONE time
  bool wasExecuted;

};

/** 
* @brief Constructor handling the new systems paramters
*
* @param ingredients_ a reference to the IngredientsType 
* @param myParameter_ a paramter to be passed
*/
template < class IngredientsType >
UpdaterExample<IngredientsType>::UpdaterExample(IngredientsType& ingredients_, uint32_t myParameter_):
    ingredients(ingredients_), myParameter(myParameter_), wasExecuted(false)
{}

/**
* @brief The initialize function handles the new systems information before execution.
*
*/
template < class IngredientsType >
void UpdaterExample<IngredientsType>::initialize(){
  std::cout << "initialize UpdaterExample" << std::endl;
  
  std::cout << "you passed the paramter "<< myParameter <<std::endl;
  
  // by default, initialise calls execute
  execute();
}

/**
* @brief Execution of your task.
*
*/
template < class IngredientsType >
bool UpdaterExample<IngredientsType>::execute(){
  if(wasExecuted)
    return true;
  
  std::cout << "execute UpdaterExample" << std::endl;

  wasExecuted=true;
  return true;
}

/**
* @brief Standard (empty) clean up.
*
*/
template < class IngredientsType >
void UpdaterExample<IngredientsType>::cleanup(){
  
}

#endif /* LEMONADE_UPDATER_YVESEXAMPLE_H */
