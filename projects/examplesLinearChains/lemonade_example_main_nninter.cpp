#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>

#include <LeMonADE/updater/UpdaterAddLinearChains.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>

int main(int argc, char* argv[])
{
  int nChains(1),chainLength(64),type1(1),type2(3),nMCS(100),nRuns(10), nSolvent((0.5*32*32*32)/8-chainLength);
  
  typedef LOKI_TYPELIST_4(
    FeatureMoleculesIO, 
    FeatureAttributes<>,
    FeatureExcludedVolumeSc<FeatureLatticePowerOfTwo<uint8_t> >,
    FeatureNNInteractionSc<FeatureLatticePowerOfTwo>) Features;
  const uint max_bonds=4;
  typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
  typedef Ingredients<Config> IngredientsType;
  IngredientsType ingredients;
    
  RandomNumberGenerators rng;
  rng.seedAll();
  
  ingredients.setBoxX(32);
  ingredients.setBoxY(32);
  ingredients.setBoxZ(32);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);

  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.setNNInteraction(1,3,0.8);
  ingredients.synchronize();
  
  TaskManager taskManager;
  taskManager.addUpdater(new UpdaterAddLinearChains<IngredientsType>(ingredients, nChains,chainLength,type1,type1),0); 
  taskManager.addUpdater(new UpdaterAddLinearChains<IngredientsType>(ingredients, nSolvent,1,type2,type2),0);
  taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_nn.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();
  
  return 0;
}