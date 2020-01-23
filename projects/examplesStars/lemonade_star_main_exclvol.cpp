#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
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

#include <LeMonADE/updater/UpdaterSimpleConnection.h>
#include "updater/UpdaterAddStars.h"
#include "updater/UpdaterSetReactivityStarEnds.h"

int main(int argc, char* argv[])
{
  int nChains(128),nChainLength(65),nBranches(4),nBranchLength(16),type1(1),nMCS(100),nRuns(10000);
  
  typedef LOKI_TYPELIST_5(
    FeatureMoleculesIO, 
    FeatureAttributes<>,
    FeatureSystemInformationStar,
    FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo < > >,
    FeatureConnectionSc
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

  ingredients.setBoxX(64);
  ingredients.setBoxY(64);
  ingredients.setBoxZ(64);
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
  taskManager.addUpdater(new UpdaterAddStars<IngredientsType>(ingredients, nChains,nChainLength,
                                                              nBranches,nBranchLength,
                                                              type1,type1),0); 

  taskManager.addUpdater(new UpdaterSetReactivityStarEnds<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,0),0); 

  taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
  taskManager.addUpdater(new UpdaterSimpleConnection<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
   
  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_ev_65_128_jan2020.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));

  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_last_65_128_jan2020.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE));

  taskManager.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients,"Rg2_65_128_jan2020.dat"),nMCS);
  taskManager.addAnalyzer(new AnalyzerBondLength<IngredientsType>(ingredients, "bondLength_65_128_jan2020.dat"),nMCS); 
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();
  
  return 0;
} 
