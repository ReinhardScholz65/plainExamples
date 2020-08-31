#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
//#include <LeMonADE/feature/FeatureConnectionScAB.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/feature/FeatureSystemInformationStars.h>

// #include <LeMonADE/updater/UpdaterAddLinearChains.h>
// #include <LeMonADE/updater/UpdaterAddStars.h>

#include <LeMonADE/updater/UpdaterSimpleSimulator.h>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include "analyzer/AnalyzerRadiusOfGyration.h"
#include "analyzer/AnalyzerBondLength.h"

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/utility/DepthIterator.h>

#include <LeMonADE/updater/UpdaterSimpleConnectionAB.h>
/7#include <LeMonADE/updater/UpdaterSimpleConnection.h>

#include "updater/UpdaterAddStarsA12B45.h"
//#include "updater/UpdaterSetReactivityStarEnds.h"
//#include "updater/UpdaterSetReactivityStarEndsAB.h"
#include "updater/UpdaterSetReactivityStarEndsA12B45.h"

int main(int argc, char* argv[])
{
  int nChainsA(100),nChainsB(100),nChainLengthA(73),nChainLengthB(97),nBranchLengthA(18),nBranchLengthB(24),nBranches(4),type1(1),type2(2),nMCS(100),nRuns(1);
  
  typedef LOKI_TYPELIST_5(
    FeatureMoleculesIO, 
    FeatureAttributes<>,
    FeatureSystemInformationStar,
    //FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo < > >,
    //FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <uint8_t> >,  /* not needed, but does not harm when using LOKI_TYPELIST_6
    //FeatureNNInteractionSc< FeatureLatticePowerOfTwo < > >,
    FeatureNNInteractionSc<FeatureLatticePowerOfTwo>,
    FeatureConnectionSc
    //FeatureConnectionScAB
    ) Features;
  const uint max_bonds=4;
  typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
  typedef Ingredients<Config> IngredientsType;
  IngredientsType ingredients;

  
  //typedef Move<Config> MoveType;
  //typedef Move MoveType;
  //MoveType move;
  
  RandomNumberGenerators rng;
  rng.seedAll();

  ingredients.setBoxX(128);
  ingredients.setBoxY(128);
  ingredients.setBoxZ(128);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();

  //ingredients.setNumStarsA(nChainsA);
  //ingredients.setNumStarsB(nChainsB);
  //ingredients.setNumChainLengthA(nChainLengthA);
  //ingredients.setNumChainLengthB(nChainLengthB);
  //ingredients.setNumBranchLengthA(nBranchLengthA);
  //ingredients.setNumBranchLengthB(nBranchLengthB);

  ingredients.setNumBranches(nBranches);

  // ingredients.setNNInteraction(1,2,0.0);
 
  ingredients.synchronize();
  
  //
  // run intialization, and run UpdaterSimpleConnection over nMCS=100, nRuns=100 (up to 10000)
  //
  
  TaskManager taskManager;

  //
  // create stars of type A, assign reactivity of branch ends
  //
  
  ingredients.setNumStars(nChainsA);
  ingredients.setNumChainLength(nChainLengthA);
  ingredients.setNumBranchLength(nBranchLengthA);

  ingredients.synchronize();

  type1=1;
  type2=2;
  taskManager.addUpdater(new UpdaterAddStarsAB<IngredientsType>(ingredients, nChainsA,nChainLengthA,nBranchLengthA,
                                                              nBranches,
                                                              type1,type2),0); 

  //taskManager.addUpdater(new UpdaterSetReactivityStarEndsAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,0),0); 

  //
  // create stars of type B, assign reactivity of branch ends
  //
  
  ingredients.setNumStars(nChainsB);
  ingredients.setNumChainLength(nChainLengthB);
  ingredients.setNumBranchLength(nBranchLengthB);
  
  ingredients.synchronize();

  type1=4;
  type2=5;
  taskManager.addUpdater(new UpdaterAddStarsAB<IngredientsType>(ingredients, nChainsB,nChainLengthB,nBranchLengthB,
                                                              nBranches,
                                                              type1,type2),0); 

  //taskManager.addUpdater(new UpdaterSetReactivityStarEndsAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,0),0); 
  
  //
  // check number of monomers
  //
  
  int nMonomers;
  nMonomers=nChainsA*nChainLengthA+nChainsB*nChainLengthB;
  std::cout << "number of monomers (A stars and B stars) " << nMonomers<< std::endl;

  //taskManager.addUpdater(new UpdaterSetReactivityStarEndsAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,0),0); 

  //taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
  ////taskManager.addUpdater(new UpdaterSimpleConnection<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  //taskManager.addUpdater(new UpdaterSimpleConnectionAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  
  taskManager.addAnalyzer(new 
  AnalyzerWriteBfmFile<IngredientsType>("config_ev_100_stars_73_97.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));

  taskManager.addAnalyzer(new
  AnalyzerWriteBfmFile<IngredientsType>("config_last_100_stars_73_97.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE));

  //taskManager.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients,"Rg2_10_stars_73_97.dat"),nMCS);
  //taskManager.addAnalyzer(new AnalyzerBondLength<IngredientsType>(ingredients, "bondLength_10_stars_73_97.dat"),nMCS); 
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();

  //
  // assign reactivity for end groups of star branches
  //
  nMCS=100;
  nRuns=1;

  TaskManager taskManager1;
  
  //! Reactivity (especially for branch ends)
  bool reactivity;
  reactivity=false;

  std::cout << "  " << std::endl;
  std::cout << "#################################  " << std::endl;
  std::cout << "  " << std::endl;
  std::cout << "system size " << ingredients.getMolecules().size() << std::endl;
  std::cout << "system size " << nMonomers << std::endl;
  std::cout << "  " << std::endl;
  std::cout << "#################################  " << std::endl;
  std::cout << "  " << std::endl;
  std::cout << "star size A " << nChainLengthA << std::endl;
  std::cout << "branch length A " << nBranchLengthA << std::endl;
  std::cout << "  " << std::endl;

  //
  // For icontrol=0, the code below identifies branch ends by careful counting through star indices and branches.
  // Then, the branch ends found are made reactive. 
  //
  // For icontrol=1, the code below identifies branch ends by checking whether the numerical number of existing bonds is 1.  
  // Then, the branch ends found are made reactive. 
  //
  // Both variants have been tested and work properly.
  //
  
  uint32_t icontrol=1;

  if(icontrol == 0) {
    
  for(uint32_t iMonomer=0;iMonomer<(nChainsA*nChainLengthA);iMonomer++){    // loop over A stars
    //std::cout << "  " << std::endl;
    //std::cout << "monomer index " << iMonomer << std::endl;
    //if(iMonomer<nChainsA*nChainLengthA) {     // loop over A stars
	uint32_t starindexA=iMonomer/nChainLengthA;
	uint32_t iMonomerInStar=iMonomer-starindexA*nChainLengthA;
	//std::cout << "star A " << starindexA << " monomer in star " << iMonomerInStar << std::endl;
	  if(nBranchLengthA*(iMonomerInStar/nBranchLengthA) == iMonomerInStar) {
	    if(iMonomerInStar == 0) { 
              std::cout << "  " << std::endl;
              std::cout << "#################################  " << std::endl;
              std::cout << "  " << std::endl;
              std::cout << "monomer index " << iMonomer << std::endl;
              std::cout << "star A " << starindexA << " star center with monomer index " << iMonomer << std::endl;
              std::cout << "number of links " <<  ingredients.getMolecules().getNumLinks(iMonomer) << std::endl;
	      if( ingredients.getMolecules()[iMonomer].isReactive() ) {
                std::cout << "reactive site, reactivity " <<  reactivity << std::endl;
	        }
	    }
	    if(iMonomerInStar != 0) { 
              std::cout << "  " << std::endl;
              std::cout << "monomer index " << iMonomer << std::endl;
              std::cout << "star A " << starindexA << " branch end with monomer index " << iMonomer << std::endl;
              std::cout << "number of links " <<  ingredients.getMolecules().getNumLinks(iMonomer) << std::endl;
              reactivity=true;
              uint32_t numMaxLinks=2;          
              ingredients.modifyMolecules()[iMonomer].setReactive(reactivity);
              ingredients.modifyMolecules()[iMonomer].setNumMaxLinks(numMaxLinks);	      
              std::cout << "maximum number of links " <<  ingredients.getMolecules()[iMonomer].getNumMaxLinks() << std::endl;
	      if( ingredients.getMolecules()[iMonomer].isReactive() ) {
                std::cout << "reactive site, reactivity " <<  reactivity << std::endl;
	        }
              //std::cout << "reactivity " <<  ingredients.getMolecules()[iMonomer].getReactive() << std::endl;
	    }
          }
	}

  } // end if(icontrol == 0)
  
  std::cout << "  " << std::endl;
  std::cout << "#################################  " << std::endl;
  std::cout << "  " << std::endl;
  std::cout << "star size B " << nChainLengthB << std::endl;
  std::cout << "branch length B " << nBranchLengthB << std::endl;

  std::cout << "  " << std::endl;
  std::cout << "#################################  " << std::endl;
  
  if(icontrol == 0) {
      
  for(uint32_t iMonomer=nChainsA*nChainLengthA;iMonomer<(nChainsA*nChainLengthA+nChainsB*nChainLengthB);iMonomer++){    // loop over B stars
    //std::cout << "  " << std::endl;
    //std::cout << "monomer index " << iMonomer << std::endl;
    //if(iMonomer<nChainsA*nChainLengthA) {     // loop over A stars
    uint32_t starindexB=(iMonomer-nChainsA*nChainLengthA)/nChainLengthB;
	uint32_t iMonomerInStar=iMonomer-nChainsA*nChainLengthA-starindexB*nChainLengthB;
	//std::cout << "star A " << starindexA << " monomer in star " << iMonomerInStar << std::endl;
	  if(nBranchLengthB*(iMonomerInStar/nBranchLengthB) == iMonomerInStar) {
	    if(iMonomerInStar == 0) { 
              std::cout << "  " << std::endl;
              std::cout << "#################################  " << std::endl;
              std::cout << "  " << std::endl;
              std::cout << "monomer index " << iMonomer << std::endl;
              std::cout << "star B " << starindexB << " star center with monomer index " << iMonomer << std::endl;
              std::cout << "number of links " <<  ingredients.getMolecules().getNumLinks(iMonomer) << std::endl;
	      if( ingredients.getMolecules()[iMonomer].isReactive() ) {
                std::cout << "reactive site, reactivity " <<  reactivity << std::endl;
	        }
	    }
	    if(iMonomerInStar != 0) { 
              std::cout << "  " << std::endl;
              std::cout << "monomer index " << iMonomer << std::endl;
              std::cout << "star B " << starindexB << " branch end with monomer index " << iMonomer << std::endl;
              std::cout << "number of links " <<  ingredients.getMolecules().getNumLinks(iMonomer) << std::endl;
	      if(ingredients.getMolecules().getNumLinks(iMonomer) == 1) {
                reactivity=true;
                uint32_t numMaxLinks=2;          
                ingredients.modifyMolecules()[iMonomer].setReactive(reactivity);
                ingredients.modifyMolecules()[iMonomer].setNumMaxLinks(numMaxLinks);
	        }
              std::cout << "maximum number of links " <<  ingredients.getMolecules()[iMonomer].getNumMaxLinks() << std::endl;
	      if( ingredients.getMolecules()[iMonomer].isReactive() ) {
                std::cout << "reactive site, reactivity " <<  reactivity << std::endl;
	        }
              //std::cout << "reactivity " <<  ingredients.getMolecules()[iMonomer].getReactive() << std::endl;
	    }
          }
	}

  } // end if(icontrol == 0)
  
  if(icontrol == 1) {

    for(uint32_t iMonomer=0;iMonomer<ingredients.getMolecules().size();iMonomer++){    // loop over all monomers
	      if(ingredients.getMolecules().getNumLinks(iMonomer) == 1) {
                reactivity=true;
                uint32_t numMaxLinks=2;          
                ingredients.modifyMolecules()[iMonomer].setReactive(reactivity);
                ingredients.modifyMolecules()[iMonomer].setNumMaxLinks(numMaxLinks);
                std::cout << "  " << std::endl;
                std::cout << "monomer index " << iMonomer << std::endl;
                std::cout << "number of links " <<  ingredients.getMolecules().getNumLinks(iMonomer) << std::endl;
                std::cout << "maximum number of links " <<  ingredients.getMolecules()[iMonomer].getNumMaxLinks() << std::endl;
                std::cout << "reactive site, reactivity " <<  reactivity << std::endl;
	        }
    }
	  
  } // end if(icontrol == 1)
  
  taskManager1.addAnalyzer(new 
  AnalyzerWriteBfmFile<IngredientsType>("config_ev_100_stars_73_97.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));

  taskManager1.addAnalyzer(new
  AnalyzerWriteBfmFile<IngredientsType>("config_last_100_stars_73_97.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE));
  
  taskManager1.initialize();
  taskManager1.run(nRuns);
  taskManager1.cleanup();
  
  return 0;

  //
  // run UpdaterSimpleConnection over nMCS=1000, nRuns=90 (10000 to 100000)
  //
  nMCS=1000;
  nRuns=50;
 
  TaskManager taskManager2;
  //taskManager2.addUpdater(new UpdaterSetReactivityStarEndsAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,0),0); 

  taskManager2.addUpdater(new UpdaterSimpleConnection<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  //taskManager.addUpdater(new UpdaterSimpleConnectionAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  
  taskManager2.addAnalyzer(new 
  AnalyzerWriteBfmFile<IngredientsType>("config_ev_100_stars_73_97.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));

  
  taskManager2.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_last_100_stars_73_97.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE));

  //taskManager2.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients,"Rg2_65_1024.dat"),nMCS);
  //taskManager2.addAnalyzer(new AnalyzerBondLength<IngredientsType>(ingredients, "bondLength_65_1024.dat"),nMCS); 
  
  taskManager2.initialize();
  taskManager2.run(nRuns);
  taskManager2.cleanup();

  return 0;
  
  //
  // run UpdaterSimpleConnection over nMCS=10000, nRuns=90 (100000 to 1000000)
  //
  nMCS=10000;
  nRuns=90;
 
  TaskManager taskManager3;
  taskManager3.addUpdater(new UpdaterSimpleConnection<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  //taskManager3.addUpdater(new UpdaterSimpleConnectionAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  
  taskManager3.addAnalyzer(new 
  AnalyzerWriteBfmFile<IngredientsType>("config_ev_65_1024.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));

  taskManager3.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_last_65_1024.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE));

  taskManager3.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients,"Rg2_65_1024.dat"),nMCS);
  taskManager3.addAnalyzer(new AnalyzerBondLength<IngredientsType>(ingredients, "bondLength_65_1024.dat"),nMCS); 
  
  taskManager3.initialize();
  taskManager3.run(nRuns);
  taskManager3.cleanup();
  
  return 0;
} 
