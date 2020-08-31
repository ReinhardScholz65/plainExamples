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

#ifndef LEMONADE_UPDATER_SETUP_CHAINSAB
#define LEMONADE_UPDATER_SETUP_CHAINSAB
/**
 * @file
 *
 * @class UpdaterAddChainsAB
 *
 * @brief Updater to create a solution of monodisperse branched stars.
 *
 * @details This is a simple implementation of a system setup starting from an empty ingredients
 * or a system with some monomers inside. This updater requires FeatureAttributes.
 * Two tags are added to the monomers in alternating manner, usually needed for GPU computing.  // not realized in this first attempt
 *                                                                                              // RS, 05 Nov 2019
 *
 * @tparam IngredientsType
 *
 * @param ingredients_ The system, holding either an empty simulation box for system setup
 * or a prefilled ingredients where the stars shall be added
 * @param NStar_ number of stars that are added to ingredients
 * @param NMonoPerStar_ number of monomers in each star
 * @param NBranchPerStar_ number of branches in each star
 * @param NMonoPerBranch_ number of monomers in each branch (excluding central monomer)
 * @param type1_ attribute tag of "even" monomers
 * @param type2_ attribute tag of "odd" monomers
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>
#include <cmath>

template<class IngredientsType>
class UpdaterAddChainsAB: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;

public:
  UpdaterAddChainsAB(IngredientsType& ingredients_, uint32_t NStar_, uint32_t NMonoPerStar_, 
                  uint32_t NBranchPerStar_, uint32_t NMonoPerBranch_, 
                  int32_t type1_=1, int32_t type2_=2, bool IsSolvent=false);

  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();

  //! getter function for write compressed solvent bool
  const bool getIsSolvent() const {return IsSolvent;}

  //! getter function for number of stars
  const int32_t getNStar() const {return NStar;}
  
  //! getter function for number of monomers in stars
  const int32_t getNMonoPerStar() const {return NMonoPerStar;}

   //! getter function for number of branches in stars
  const int32_t getNBranchPerStar() const {return NBranchPerStar;}

  //! getter function for number of monomers in branch
  const int32_t getNMonoPerBranch() const {return NMonoPerBranch;}

  //! getter function for calculated density
  const double getDensity() const {return density;}

private:
  // provide access to functions of UpdaterAbstractCreate used in this updater
  using BaseClass::ingredients;
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::linearizeSystem;

  //! number of stars in the box
  uint32_t NStar;
  
  //! number of monomers in a star
  uint32_t NMonoPerStar;

  //! number of branches in a star
  uint32_t NBranchPerStar;

  //! number of monomers in a branch
  uint32_t NMonoPerBranch;  

  //! lattice occupation density
  double density;

  //! bool for execution
  bool wasExecuted;

  //! attribute tag of even monomers
  int32_t type1;

  //! getAttributeTag of odd monomers
  int32_t type2;

  //! bool to check if chains of size 1 should be compressed to solvent
  bool IsSolvent;
};

/**
* @brief Constructor handling the new systems paramters
* 
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param NStar_ number of stars to be added in the system instead of solvent
* @param NMonoPerStar_ number of monomers in each star
* @param NBranchPerStar_ number of branches in each star
* @param NMonoPerBranch_ number of monomers in each branch (excluding central monomer)
*/


template < class IngredientsType >
UpdaterAddChainsAB<IngredientsType>::UpdaterAddChainsAB(IngredientsType& ingredients_, uint32_t NStar_, 
                                                  uint32_t NMonoPerStar_, uint32_t NBranchPerStar_, uint32_t NMonoPerBranch_, 
                                                  int32_t type1_, int32_t type2_, bool IsSolvent_):
BaseClass(ingredients_), NStar(NStar_), 
NMonoPerStar(NMonoPerStar_), NBranchPerStar(NBranchPerStar_), NMonoPerBranch(NMonoPerBranch_), 
density(0.0), wasExecuted(false), type1(type1_), type2(type2_), IsSolvent(IsSolvent_)
{}

/**
* @brief initialise function, calculate the target density to compare with at the end.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddChainsAB<IngredientsType>::initialize(){
  std::cout << "initialize UpdaterAddChainsAB" << std::endl;

  // get the target density from the sum of existing monomers and the new added chains
  density=(double)( ingredients.getMolecules().size() + NMonoPerStar*NStar ) * 8  /(double)( ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ() );

  std::cout << "add "<<NStar*NMonoPerStar<<" monomers to the box"<<std::endl;

  execute();
}

/**
* @brief Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterAddChainsAB<IngredientsType>::execute(){
  if(wasExecuted)
    return true;

  std::cout << "execute UpdaterAddChainsAB" << std::endl;

  //loop over stars and star monomers and build it up
  for(uint32_t i=0;i<(NStar);i++){
    for(uint32_t j=0;j<(NMonoPerStar);j++){
      //
      if(j<(NMonoPerStar/2)){        // first half of chain: monomers of type A, type1=1 or type1=2
        if(2*(j/2)==j){              // monomer with even index 
          type1 = 2;                 // A, itype=2      
	  if(j==0){
            addSingleMonomer(type1); // first A monomer
	  }   
          else{
            addMonomerToParent(ingredients.getMolecules().size()-1,type1);  // size()-1 shall correspond to previous monomer
	  }
	}
        if(2*(j/2)!=j){              // monomer with odd index       
          type1 = 1;                 // A, itype=1
          addMonomerToParent(ingredients.getMolecules().size()-1,type1);  // size()-1 shall correspond to previous monomer
        }
	}  // end first half chain
        //
      if(j>=(NMonoPerStar/2)){       // second half of chain: monomers of type B, type1=4 or type1=5
        if(2*(j/2)==j){              // monomer with even index       
          type1 = 5;                 // B, itype=5
          addMonomerToParent(ingredients.getMolecules().size()-1,type1);  // size()-1 shall correspond to previous monomer
	}
        if(2*(j/2)!=j){              // monomer with odd index       
          type1 = 4;                 // B, itype=4
          addMonomerToParent(ingredients.getMolecules().size()-1,type1);  // size()-1 shall correspond to previous monomer
        }
	}  // end second half of chain
	//
    } // end loop j
  }  // end loop i
	
    

  ingredients.synchronize();
  double lattice_volume(ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());
  if(std::abs(density - ( (double)(ingredients.getMolecules().size()*8) / lattice_volume )) > 0.0000000001 ){
    std::cout << density << " " <<( (ingredients.getMolecules().size()*8) / lattice_volume)<<std::endl;
    throw std::runtime_error("UpdaterAddChainsAB: number of monomers in molecules does not match the calculated number of monomers!");
  }else{
    std::cout << "real lattice occupation density =" << (8*ingredients.getMolecules().size()) / lattice_volume<<std::endl;
    wasExecuted=true;
    // if we allow for solvent compression AND added single monomers (solvent): compress solvent
    if(IsSolvent && NMonoPerStar==1){
      ingredients.setCompressedOutputIndices(ingredients.getMolecules().size()-(NMonoPerStar*NStar),ingredients.getMolecules().size()-1);
    }else{
      linearizeSystem();
    }
    return true;
  }
}

/**
* @brief Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddChainsAB<IngredientsType>::cleanup(){

}


#endif /* LEMONADE_UPDATER_SETUP_CHAINSAB */
