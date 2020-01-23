#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
//#include <LeMonADE/feature/FeatureConnectionScAB.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
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

//#include <LeMonADE/updater/UpdaterSimpleConnectionAB.h>
#include <LeMonADE/updater/UpdaterSimpleConnection.h>

#include "updater/UpdaterAddStarsAB.h"
#include "updater/UpdaterSetReactivityStarEnds.h"
//#include "updater/UpdaterSetReactivityStarEndsAB.h"

int main(int argc, char* argv[])
{
  int nChains(1024),nChainLength(65),nBranches(4),nBranchLength(16),type1(1),type2(2),nMCS(100),nRuns(100000);
  
  typedef LOKI_TYPELIST_5(
    FeatureMoleculesIO, 
    FeatureAttributes<>,
    FeatureSystemInformationStar,
    FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo < > >,
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
  ingredients.setNumStars(nChains);
  ingredients.setNumBranches(nBranches);
  ingredients.setNumChainLength(nChainLength);
  ingredients.setNumBranchLength(nBranchLength);
  ingredients.synchronize();
  
  TaskManager taskManager;
  taskManager.addUpdater(new UpdaterAddStarsAB<IngredientsType>(ingredients, nChains,nChainLength,
                                                              nBranches,nBranchLength,
                                                              type1,type2),0); 

  taskManager.addUpdater(new UpdaterSetReactivityStarEnds<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,0),0); 
  //taskManager.addUpdater(new UpdaterSetReactivityStarEndsAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,0),0); 

  //taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
  taskManager.addUpdater(new UpdaterSimpleConnection<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  //taskManager.addUpdater(new UpdaterSimpleConnectionAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  
  //taskManager.addAnalyzer(new 
  //AnalyzerWriteBfmFile<IngredientsType>("config_ev_65_1024.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));

  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_last_65_1024.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE));

  taskManager.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients,"Rg2_65_1024.dat"),nMCS);
  taskManager.addAnalyzer(new AnalyzerBondLength<IngredientsType>(ingredients, "bondLength_65_1024.dat"),nMCS); 
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();
  
  return 0;
} 
