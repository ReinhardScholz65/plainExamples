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

#ifndef LEMONADE_UPDATER_SETREACTIVITYSTARENDS_H
#define LEMONADE_UPDATER_SETREACTIVITYSTARENDS_H


#include <LeMonADE/updater/moves/MoveLocalBase.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>

/**
 * @file
 *
 * @class UpdaterSetReactivityStarEnds
 *
 * @brief Simple simulation updater for generic polyrotaxanes.
 *
 * @details It takes the type of move as template argument MoveType
 * and the number of mcs to be executed as argument for the constructor
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 * @tparam MoveType name of the specialized move.
 */

template<class IngredientsType,class MoveType, class ConnectionMoveType>
class UpdaterSetReactivityStarEnds:public AbstractUpdater
{

public:
  /**
   * @brief Standard Constructor initialized with ref to Ingredients and MCS per cycle
   *
   * @param ing a reference to the IngredientsType - mainly the system
   * @param steps MCS per cycle to performed by execute()
   */
  
  /**
   * additonal parameters needed to administrate number of stars, number of branches, etc.:
   * 
   * @param NStar_ number of stars in ingredients
   * @param NMonoPerStar_ number of monomers in each star
   * @param NBranchPerStar_ number of branches in each star
   * @param NMonoPerBranch_ number of monomers in each branch (excluding central monomer)
   * @param type1_ attribute tag of "even" monomers
   * @param type2_ attribute tag of "odd" monomers
   **/

//  UpdaterSetReactivityStarEnds(IngredientsType& ing,uint32_t steps = 1 )
  UpdaterSetReactivityStarEnds(IngredientsType& ing,uint32_t steps = 1 )
  :ingredients(ing),nsteps(steps),NReactedSites(0),NReactiveSites(0){}


//UpdaterSetReactivityStarEnds(IngredientsType& ing,uint32_t steps = 1,uint32_t NStar_, uint32_t NMonoPerStar_, 
//                    uint32_t NBranchPerStar_, uint32_t NMonoPerBranch_, 
//                    int32_t type1_=1, int32_t type2_=2, bool IsSolvent=false,
//                    uint32_t NReactedSites,uint32_t NReactiveSites )
//  :ingredients(ing),nsteps(steps),
//                    uint32_t NStar_, uint32_t NMonoPerStar_, 
//                    uint32_t NBranchPerStar_, uint32_t NMonoPerBranch_, 
//                    int32_t type1_=1, int32_t type2_=2, bool IsSolvent=false,
//                    uint32_t NReactedSites(0),uint32_t NReactiveSites(0){}

  /**
   * @brief This checks all used Feature and applies all Feature if all conditions are met.
   *
   * @details This function runs over \a steps MCS and performs the moves.
   * It setting the age of the system and prints a simple simple simulation speed
   * in the number of attempted monomer moves per s (tried and performed monomer moves/s).
   *
   * @return True if function are done.
   */
  bool execute();

  /**
   * @brief This function is called \a once in the beginning of the TaskManager.
   *
   * @details It´s a virtual function for inheritance.
   * Use this function for initializing tasks (e.g. init SDL)
   *
   **/
  virtual void initialize();

  /**
   * @brief This function is called \a once in the end of the TaskManager.
   *
   * @details It´s a virtual function for inheritance.
   * Use this function for cleaning tasks (e.g. destroying arrays, file outut)
   *
   **/
  virtual void cleanup(){};

  
  //! get conversion of the reaction process
  double getConversion(){return (double(NReactedSites))/(double(NReactiveSites));};
  
protected:
  
  //! A reference to the IngredientsType - mainly the system
  IngredientsType& ingredients;
  
  //! number of monomers which could have a reaction bond
  uint32_t NReactiveSites;

  //! number of already reacted monomers
  uint32_t NReactedSites;
  
  //! Specialized move to be used for the movement of the monomers 
  MoveType move;

  //! Specialized move to be used for the connection between reactive monomers
  ConnectionMoveType connectionMove;  
  
  //! random number generator (seed set in main program)
  RandomNumberGenerators rng;
  
private:

  //! Number of mcs to be executed
  uint32_t nsteps;

  //! Number of stars
  uint32_t NStar;

  //! Number of monomers per star
  uint32_t NMonoPerStar;

  //! Number of branches per star
  uint32_t NBranchPerStar;
 
   //! Number of monomers per branch
  uint32_t NMonoPerBranch;
 
 //! Reactivity (especially for branch ends)
  bool reactivity;
  
  //! Number of maximimum possible links/connections to another not-yet-connected monomers
  uint32_t numMaxLinks;
 
};
/**Implementation of the member functions
 * @brief 
 */
template<class IngredientsType,class MoveType, class ConnectionMoveType>
bool UpdaterSetReactivityStarEnds<IngredientsType,MoveType,ConnectionMoveType>::execute()
{
  
// 	time_t startTimer = time(NULL); //in seconds
// 	std::cout<<"connection mcs "<<ingredients.getMolecules().getAge() << " passed time " << ((difftime(time(NULL), startTimer)) ) <<std::endl;

	for(int n=0;n<nsteps;n++)
	{
		// std::vector<uint32_t> TagedMonomers;

		for(size_t m=0;m<ingredients.getMolecules().size();m++)
		{
			// if( ingredients.getMolecules()[m].isReactive() )
			//	 TagedMonomers.push_back(m);

			move.init(ingredients);

			if(move.check(ingredients)==true)
			{
				move.apply(ingredients);
			}
			else // move is reject due to e.g. excluded volume, bond length, Metropolis etc
			{
				// as we connect as face-to-face colliding algorithm
				// use the move direction in previous check as destination direction
				if( ingredients.getMolecules()[move.getIndex()].isReactive() )
				{
					// collision reaction
					connectionMove.init(ingredients,move.getIndex(), 2*move.getDir());

					if ( connectionMove.check(ingredients) )
					{
						connectionMove.apply(ingredients);
						NReactedSites+=2; // to check
					}
				}
			}
		}

		ingredients.modifyMolecules().setAge(ingredients.getMolecules().getAge()+1);
	}
// 	std::cout <<"Conversion at "<<ingredients.getMolecules().getAge() << " is " << getConversion()  << std::endl;
// 	std::cout<<"connection mcs "<<ingredients.getMolecules().getAge() << " with " << (((1.0*nsteps)*ingredients.getMolecules().size())/(difftime(time(NULL), startTimer)) ) << " [attempted connections/s]" <<std::endl;
// 	std::cout<<"connection mcs "<<ingredients.getMolecules().getAge() << " passed time " << ((difftime(time(NULL), startTimer)) ) << " with " << nsteps << "connection MCS "<<std::endl;
	return false;
};

template<class IngredientsType,class MoveType, class ConnectionMoveType>
void  UpdaterSetReactivityStarEnds<IngredientsType,MoveType,ConnectionMoveType>::initialize()
{
    // attempt to initialize reactivity of branch ends
    //

    {
        
    NStar=ingredients.getNumStars();
    NBranchPerStar=ingredients.getNumBranches();
    NMonoPerStar=ingredients.getNumChainLength();
    NMonoPerBranch=ingredients.getNumBranchLength();

 	std::cout << NStar <<" stars "<< std::endl;
  	std::cout << NBranchPerStar <<" branches per star "<< std::endl;
  	std::cout << NMonoPerStar <<" monomers per star "<< std::endl;
   	std::cout << NMonoPerBranch <<" monomers per branch "<< std::endl;
   	std::cout << std::endl;
 
    for(uint32_t istar=0;istar<(NStar);istar++){
       for(uint32_t j=0;j<(NMonoPerStar);j++){
            reactivity=false;
            ingredients.modifyMolecules()[istar*NMonoPerStar+j].setReactive(reactivity);
        	//std::cout << istar*NMonoPerStar+j <<" monomer "<< std::endl;
       }
    }
    
    for(uint32_t istar=0;istar<(NStar);istar++){
        for(uint32_t ibranch=0;ibranch<(NBranchPerStar);ibranch++){
            reactivity=true;
            numMaxLinks=2;
            ingredients.modifyMolecules()[istar*NMonoPerStar+(ibranch+1)*NMonoPerBranch].setReactive(reactivity);
            ingredients.modifyMolecules()[istar*NMonoPerStar+(ibranch+1)*NMonoPerBranch].setNumMaxLinks(numMaxLinks);       
        	//std::cout << istar <<" star number; " << istar*NMonoPerStar+(ibranch+1)*NMonoPerBranch <<" monomer number "<<std::endl;            
        }
    }

    //for(uint32_t istar=0;istar<(NStar);istar++){
    //   for(uint32_t j=0;j<(NMonoPerStar);j++){
    //        reactivity=false;
    //        ingredients.modifyMolecules()[istar*NMonoPerStar+j].setReactive(reactivity);
    //        if(j>=NMonoPerStar-NBranchPerStar){
    //            reactivity=true;
    //            numMaxLinks=2;
    //            ingredients.modifyMolecules()[istar*NMonoPerStar+j].setReactive(reactivity);
    //            ingredients.modifyMolecules()[istar*NMonoPerStar+j].setNumMaxLinks(numMaxLinks);
    //        }
    //    }
    // }
    
    }
    
    //
    //end initialize reactivity

	for(size_t i = 0 ; i < ingredients.getMolecules().size(); i++ )
	{
		if ( ingredients.getMolecules()[i].isReactive() )
		{
			uint32_t NLinks(ingredients.getMolecules().getNumLinks(i));
			uint32_t nIrreversibleBonds=0;
			for (uint32_t n = 0 ; n < NLinks ;n++)
			{
				uint32_t neighbor(ingredients.getMolecules().getNeighborIdx(i,n));
				if( ingredients.getMolecules()[neighbor].isReactive() )
					NReactedSites++;
				else
					nIrreversibleBonds++;
			}
			NReactiveSites+=(ingredients.getMolecules()[i].getNumMaxLinks()-nIrreversibleBonds);
		}
	}
    ingredients.synchronize();
	std::cout << NReactedSites <<" "<<NReactiveSites<<" "<<getConversion()<<std::endl;

	if ( NReactiveSites == 0  )
	{
	  std::stringstream errormessage;
	  errormessage << "UpdaterSimpleConnection::initialize(): The number of possible reactive sites is zero. \n"
		       << "Check if the system is correctly setup with reactivity!";
	  throw std::runtime_error(errormessage.str());
	} 
};
#endif 	/*LEMONADE_UPDATER_SETREACTIVITYSTARENDS_H*/

