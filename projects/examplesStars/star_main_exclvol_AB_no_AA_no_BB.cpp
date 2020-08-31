//
// intended: AA and BB reactions suppressed
//
// realized:
//
// with this combination of include commands, the type of star alternates (1-65: A, 66-130: B, ...)
//
// bonds can only be formed between stars of different type A and B
//
// AB and BA bonds (i.e. between different stars) are allowed
//
// AA and BB bonds (i.e. between two A stars, between two B stars, or within a star of type A or type B) cannot be formed
//
// most important line:
//
// last variable of updater has to be the source code suppressing the possibility of AA or BB bonds
// UpdaterSimpleConnectionAB<IngredientsType,MoveLocalSc,MoveConnectScAB>
// where MoveConnectScAB.h inhibits these AA or BB bonds by setting 'false' when tagAttributes of chain ends coincide 
//

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

#include <LeMonADE/updater/UpdaterSimpleConnectionAB.h>
//#include <LeMonADE/updater/UpdaterSimpleConnection.h>

#include "updater/UpdaterAddStarsAB.h"
//#include "updater/UpdaterSetReactivityStarEnds.h"
#include "updater/UpdaterSetReactivityStarEndsAB.h"

int main(int argc, char* argv[])
{
  int nChains(1024),nChainLength(65),nBranches(4),nBranchLength(16),type1(1),type2(2),nMCS(100),nRuns(100);
  
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
  
  //
  // run intialization, and run UpdaterSimpleConnection over nMCS=100, nRuns=100 (up to 10000)
  //
  
  TaskManager taskManager;
  taskManager.addUpdater(new UpdaterAddStarsAB<IngredientsType>(ingredients, nChains,nChainLength,
                                                              nBranches,nBranchLength,
                                                              type1,type2),0); 

  //taskManager.addUpdater(new UpdaterSetReactivityStarEnds<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,0),0); 
  taskManager.addUpdater(new UpdaterSetReactivityStarEndsAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,0),0); 

  //taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
  //taskManager.addUpdater(new UpdaterSimpleConnection<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  taskManager.addUpdater(new UpdaterSimpleConnectionAB<IngredientsType,MoveLocalSc,MoveConnectScAB>(ingredients,nMCS));
  
  taskManager.addAnalyzer(new 
  AnalyzerWriteBfmFile<IngredientsType>("config_ev_65_1024.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));

  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_last_65_1024.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE));

  taskManager.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients,"Rg2_65_1024.dat"),nMCS);
  taskManager.addAnalyzer(new AnalyzerBondLength<IngredientsType>(ingredients, "bondLength_65_1024.dat"),nMCS); 
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();

  return 0;
   
  //
  // run UpdaterSimpleConnection over nMCS=1000, nRuns=90 (10000 to 100000)
  //
  nMCS=1000;
  nRuns=90;
 
  TaskManager taskManager2;
  //taskManager2.addUpdater(new UpdaterSimpleConnection<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  taskManager2.addUpdater(new UpdaterSimpleConnectionAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  
  taskManager2.addAnalyzer(new 
  AnalyzerWriteBfmFile<IngredientsType>("config_ev_65_1024.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));

  taskManager2.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_last_65_1024.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE));

  taskManager2.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients,"Rg2_65_1024.dat"),nMCS);
  taskManager2.addAnalyzer(new AnalyzerBondLength<IngredientsType>(ingredients, "bondLength_65_1024.dat"),nMCS); 
  
  taskManager2.initialize();
  taskManager2.run(nRuns);
  taskManager2.cleanup();
  
  //
  // run UpdaterSimpleConnection over nMCS=10000, nRuns=90 (100000 to 1000000)
  //
  nMCS=10000;
  nRuns=90;
 
  TaskManager taskManager3;
  //taskManager3.addUpdater(new UpdaterSimpleConnection<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  taskManager3.addUpdater(new UpdaterSimpleConnectionAB<IngredientsType,MoveLocalSc,MoveConnectSc>(ingredients,nMCS));
  
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
